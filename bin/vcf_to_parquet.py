from cyvcf2 import VCF, Variant
import polars as pl
from typing import Dict, List, Any, Optional, Tuple
import numpy as np
import argparse
import os


def vcf_to_parquet(vcf_path: str, parquet_path: str, output_metadata_schema = False, include_sample_col = True) -> None:
    """
    Convert a single-sample VCF file to Parquet format.
    
    Parameters
    ----------
    vcf_path : str
        Path to input VCF or VCF.gz file
    parquet_path : str
        Path to output Parquet file
    output_metadata_schema : bool, optional
        If True, output a TSV file with metadata schema for INFO and FORMAT fields
    include_sample_col : bool, optional
        If True, include a 'sample' column with the sample name from the VCF
    Notes
    -----
    - Multi-allelic sites are kept as single rows with ALT as List[Utf8]
    - INFO and FORMAT fields are stored as properly typed Polars Structs
    - Missing values (.) are converted to nulls
    - FILTER is stored as List[Utf8], and PASS is represented as null
    - INFO fields with Number=A/R/G/. are stored as Polars List types
    - FORMAT fields are extracted for the single sample in the VCF
    """
    vcf = VCF(vcf_path)
    
    # Validate single sample
    if len(vcf.samples) != 1:
        raise ValueError(f"Expected single sample VCF, got {len(vcf.samples)} samples")
    
    sample_name = vcf.samples[0]
    
    # Build schema from VCF header
    info_defs, info_meta = _parse_header_fields(vcf, 'INFO')
    format_defs, format_meta = _parse_header_fields(vcf, 'FORMAT')

    
    info_schema = _build_struct_schema(info_defs)
    format_schema = _build_struct_schema(format_defs)
    
    # Define full schema
    schema = {
        'CHROM': pl.Utf8,
        'POS': pl.Int64,
        'ID': pl.Utf8,
        'REF': pl.Utf8,
        'ALT': pl.List(pl.Utf8),
        'QUAL': pl.Float64,
        'FILTER': pl.List(pl.Utf8),
        'INFO': pl.Struct(info_schema),
        'FORMAT': pl.Struct(format_schema),
    }
    
    if include_sample_col:
        schema['sample'] = pl.Utf8
    # Collect data
    data = []
    for variant in vcf:
        row = _parse_variant(variant, info_defs, format_defs)
        if include_sample_col:
            row['sample'] = sample_name
        data.append(row)
    
    vcf.close()
    
    if output_metadata_schema:
        info_meta_df = pl.DataFrame({
            'column': list(info_meta.keys()),
            'description': list(info_meta.values()),
            'table': 'info'
        })
        format_meta_df = pl.DataFrame({
            'column': list(format_meta.keys()),
            'description': list(format_meta.values()),
            'table': 'format'
        })
        meta_df = pl.concat([info_meta_df, format_meta_df])
        meta_out_path = parquet_path.replace('.parquet', '_metadata_schema.tsv')
        meta_df.write_csv(meta_out_path, separator='\t')


    # Create DataFrame with explicit schema
    df = pl.DataFrame(data, schema=schema)
    print("Sample of VCF data parsed:")
    print(df.head())
    df.write_parquet(parquet_path, compression='snappy')
    print(f"Wrote {len(df)} variants to {parquet_path}")


def _parse_header_fields(vcf: VCF, field_type: str) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str]]:
    """Extract INFO or FORMAT field definitions from VCF header, alongside their descriptions."""
    fields = {}
    description = {}
    for row in vcf.header_iter():
        if row.type == field_type:
            field_id = row['ID']
            if field_id == 'GT':
                # cyvcf2 access genotypes as an list of form [ref, alt, phased]
                fields[field_id] = {
                    'number': '3',
                    'type': 'Integer',
                }
            else:
                fields[field_id] = {
                    'number': row['Number'],
                    'type': row['Type'],
                }
            description[field_id] = row['Description']
    return fields, description


def _build_struct_schema(field_defs: Dict[str, Dict[str, str]]) -> Dict[str, pl.DataType]:
    """Build Polars struct schema from VCF field definitions."""
    schema = {}
    for field_id, field_def in field_defs.items():
        schema[field_id] = _get_polars_type(field_def)
    return schema


def _get_polars_type(field_def: Dict[str, str]) -> pl.DataType:
    """Convert VCF type definition to Polars type."""
    vcf_type = field_def['type']
    number = field_def['number']
    
    # Determine base type
    base_type = pl.String()
    
    if vcf_type == 'Integer':
        base_type = pl.Int64()
    elif vcf_type == 'Float':
        base_type = pl.Float64()
    elif vcf_type == 'Flag':
        base_type = pl.Boolean()
    elif vcf_type == 'String':
        base_type = pl.String()
    
    # Check if should be array
    # Number can be: A (alt alleles), R (ref + alt), G (genotypes), . (variable), or integer
    if number in ['A', 'R', 'G', '.'] or (number.isdigit() and int(number) > 1):
        return pl.List(base_type)
   
    return base_type


def _parse_variant(variant: Variant, info_defs: Dict, format_defs: Dict) -> Dict[str, Any]:
    """Parse a single variant into a dictionary with proper Python types."""
    
    # Basic fields

    row = {
        'CHROM': variant.CHROM,
        'POS': variant.POS,
        'ID': variant.ID if variant.ID else None,
        'REF': variant.REF,
        'ALT': list(variant.ALT) if variant.ALT else [],
        'QUAL': float(variant.QUAL) if variant.QUAL is not None else None,
        'FILTER': _parse_filter(variant.FILTER),
    }
    

    # Parse INFO fields
    info_dict = {}
    for field_id, field_def in info_defs.items():
        value = variant.INFO.get(field_id)
        info_dict[field_id] = _convert_value(value, field_def)
    row['INFO'] = info_dict
    
    # Parse FORMAT fields
    format_dict = {}
    for field_id, field_def in format_defs.items():
        try:
            if field_id == 'GT':
                fmt_array = variant.genotypes 
                format_dict[field_id] = fmt_array[0]
            else:
                fmt_array = variant.format(field_id)
                if fmt_array is not None and len(fmt_array) > 0:
                    # Extract single sample (index 0)
                    value = fmt_array[0]
                    format_dict[field_id] = _convert_format_value(value, field_def)
                else:
                    format_dict[field_id] = None
        except (KeyError, AttributeError):
                format_dict[field_id] = None
    
    row['FORMAT'] = format_dict
    
    return row

def _parse_filter(filter_value) -> List[str] | None:
    """Parse FILTER field into list of filter strings."""
    if filter_value is None or filter_value == '.':
        return None
    # If already a list, return it
    if isinstance(filter_value, list):
        return filter_value
    
    # If string, split by semicolon (multiple filters)
    if isinstance(filter_value, str):
        if ';' in filter_value:
            return filter_value.split(';')
        else:
            # Single filter like 'PASS' or '.'
            return [filter_value] if filter_value != '.' else []
    
    # Fallback
    return [str(filter_value)]

def _convert_value(value, field_def: Dict) -> Any:
    """Convert INFO field value to Python type (not numpy)."""
    if value is None:
        return None
    
    vcf_type = field_def['type']
    number = field_def['number']
    
    # Handle Flag type
    if vcf_type == 'Flag':
        return bool(value)
    
    # Handle missing values
    if value == '.' or value == '':
        return None
    
    # Convert based on type
    if isinstance(value, (tuple, list)):
        # Array field
        converted = []
        for v in value:
            if v is None or v == '.' or v == '':
                converted.append(None)
            else:
                converted.append(_convert_scalar(v, vcf_type))
        return converted
    # when field is single value we want to maintain structure
    elif number in ['A', 'R', 'G', '.'] or (number.isdigit() and int(number) > 1):
        converted = []
        converted.append(_convert_scalar(value, vcf_type))
        return converted
    else:
        # Scalar field
        return _convert_scalar(value, vcf_type)


def _convert_format_value(value, field_def: Dict) -> Any:
    """Convert FORMAT field value from numpy array to Python types."""
    if value is None:
        return None
    
    vcf_type = field_def['type']
    number = field_def['number']
    

    # Handle numpy arrays - convert to Python types
    if isinstance(value, np.ndarray):
        # Check if it's a single element
        if value.size == 1:
            val = value.item()
            if _is_missing(val, vcf_type):
                return None
            return _convert_scalar(val, vcf_type)
        else:
            # Multiple values - convert to Python list
            converted = []
            for v in value.flat:  # Use flat to handle any shape
                if _is_missing(v, vcf_type):
                    converted.append(None)
                else:
                    converted.append(_convert_scalar(v, vcf_type))
            return converted
    
    # Handle scalar
    if _is_missing(value, vcf_type):
        return None
    return _convert_scalar(value, vcf_type)


def _is_missing(value, vcf_type: str) -> bool:
    """Check if value represents missing data."""
    if value is None:
        return True
    if isinstance(value, float) and np.isnan(value):
        return True
    if vcf_type == 'String' and value == '.':
        return True
    if value == '':
        return True
    return False


def _convert_scalar(value, vcf_type: str) -> Any:
    """Convert a single value to appropriate Python type (not numpy)."""
    if value is None or value == '.' or value == '':
        return None
    
    if vcf_type == 'Integer':
        return int(value)
    elif vcf_type == 'Float':
        return float(value)
    elif vcf_type == 'String':
        return str(value)
    elif vcf_type == 'Flag':
        return bool(value)
    else:
        return value

def main():
    parser = argparse.ArgumentParser(
        description='Convert a single-sample VCF to Parquet format'
    )
    parser.add_argument('vcf_path', help='Path to VCF file')
    parser.add_argument('--metadata-schema', action='store_true',
                        help='Output TSV file with metadata schema for INFO and FORMAT fields')
    parser.add_argument('--include-sample-col', action='store_true',
                        help='Include a "sample" column with the sample name from the VCF')
    parser.add_argument('-o', '--output', 
                        help='Output parquet file (default: basenames(vcf_path).parquet)')
    
    args = parser.parse_args()
    
    output_file = args.output if args.output else f'{os.path.basename(args.vcf_path)}.parquet'
    vcf_to_parquet(args.vcf_path, output_file, args.metadata_schema, args.include_sample_col)


if __name__ == "__main__":
    main()
