"""Command line interface for HamDiv molecular diversity calculations."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List, Optional

import click
import pandas as pd
from rdkit import Chem
from rich.console import Console
from rich.progress import track
from rich.table import Table

from .diversity import diversity_all, HamDiv, diversity_vector, diversity_vector_array

console = Console()


def print_version(ctx, param, value):
    """Print version information."""
    if not value or ctx.resilient_parsing:
        return
    
    try:
        from importlib.metadata import version
        hamdiv_version = version('hamdiv')
    except ImportError:
        hamdiv_version = "Unknown"
    
    click.echo(f"HamDiv v{hamdiv_version}")
    click.echo(f"Python: {sys.version.split()[0]}")
    click.echo(f"Platform: {sys.platform}")
    ctx.exit()


@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True, help='Show version and exit.')
@click.pass_context
def cli(ctx):
    """HamDiv - Hamiltonian diversity and molecular diversity metrics.
    
    This tool provides access to Hamiltonian diversity and all major
    molecular diversity metrics for chemical datasets.
    """
    ctx.ensure_object(dict)


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o", "--output", 
    type=click.Path(path_type=Path),
    help="Output file for diversity metrics (CSV format). If not specified, prints to stdout."
)
@click.option(
    "--format", "input_format",
    type=click.Choice(['smi', 'smiles', 'sdf'], case_sensitive=False),
    help="Input file format (auto-detected if not specified)"
)
@click.option(
    "--metric",
    type=click.Choice([
        'all', 'hamdiv', 'richness', 'intdiv', 'sumdiv', 'diam', 'sumdiam',
        'bot', 'sumbot', 'dpp', 'ncircles', 'fg', 'rs', 'bm'
    ], case_sensitive=False),
    default='all',
    help="Diversity metric to calculate (default: all)"
)
@click.option(
    "--threshold",
    type=float,
    default=0.7,
    help="Threshold for NCircles metric (default: 0.7)"
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Verbose output"
)
def calculate(
    input_file: Path,
    output: Optional[Path],
    input_format: Optional[str],
    metric: str,
    threshold: float,
    verbose: bool
):
    """Calculate molecular diversity metrics for molecules in input file.
    
    Examples:
    \b
        hamdiv calculate molecules.smi
        hamdiv calculate data.sdf --metric hamdiv --verbose
        hamdiv calculate compounds.smi --metric ncircles --threshold 0.8
        hamdiv calculate dataset.smi -o diversity_results.csv
    """
    if verbose:
        console.print(f"[green]Reading molecules from:[/green] {input_file}")
    
    # Auto-detect format
    if not input_format:
        suffix = input_file.suffix.lower()
        if suffix in ['.smi', '.smiles']:
            input_format = 'smiles'
        elif suffix == '.sdf':
            input_format = 'sdf'
        else:
            input_format = 'smiles'  # Default assumption
    
    # Read molecules
    try:
        smiles_list = read_molecules(input_file, input_format, verbose)
        if not smiles_list:
            console.print("[red]Error:[/red] No valid molecules found in input file")
            sys.exit(1)
            
        if verbose:
            console.print(f"[green]Loaded {len(smiles_list)} molecules[/green]")
        
        # Calculate diversity metrics
        if verbose:
            console.print("[yellow]Calculating diversity metrics...[/yellow]")
        
        if metric == 'all':
            # Calculate all metrics
            results = diversity_vector(smiles=smiles_list)
            if verbose:
                console.print("[green]Calculated all diversity metrics[/green]")
        else:
            # Calculate specific metric
            if metric == 'ncircles':
                metric_name = f"NCircles-{threshold}"
            else:
                metric_map = {
                    'hamdiv': 'HamDiv',
                    'richness': 'Richness',
                    'intdiv': 'IntDiv',
                    'sumdiv': 'SumDiv',
                    'diam': 'Diam',
                    'sumdiam': 'SumDiam',
                    'bot': 'Bot',
                    'sumbot': 'SumBot',
                    'dpp': 'DPP',
                    'fg': 'FG',
                    'rs': 'RS',
                    'bm': 'BM'
                }
                metric_name = metric_map[metric]
            
            value = diversity_all(smiles=smiles_list, mode=metric_name)
            results = {metric: value}
            
            if verbose:
                console.print(f"[green]Calculated {metric} = {value:.4f}[/green]")
        
        # Output results
        if output:
            if metric == 'all':
                # Save as CSV with all metrics
                df = pd.DataFrame([results])
                df.to_csv(output, index=False)
            else:
                # Save single metric
                df = pd.DataFrame([{
                    'dataset': input_file.stem,
                    'metric': metric,
                    'value': list(results.values())[0]
                }])
                df.to_csv(output, index=False)
            
            if verbose:
                console.print(f"[green]Results saved to:[/green] {output}")
        else:
            # Print to stdout
            if metric == 'all':
                # Pretty table for all metrics
                table = Table(title="Molecular Diversity Metrics")
                table.add_column("Metric", style="cyan", no_wrap=True)
                table.add_column("Value", style="white", justify="right")
                
                for metric_name, value in results.items():
                    if isinstance(value, (int, float)):
                        table.add_row(metric_name.replace('_', ' ').title(), f"{value:.4f}")
                    else:
                        table.add_row(metric_name.replace('_', ' ').title(), str(value))
                
                console.print(table)
            else:
                # Simple output for single metric
                value = list(results.values())[0]
                console.print(f"{metric}: {value:.4f}")
            
    except Exception as e:
        console.print(f"[red]Error calculating diversity:[/red] {e}")
        sys.exit(1)


@cli.command()
@click.argument("smiles_list", nargs=-1, required=True)
@click.option(
    "--metric",
    type=click.Choice([
        'all', 'hamdiv', 'richness', 'intdiv', 'sumdiv', 'diam', 'sumdiam',
        'bot', 'sumbot', 'dpp', 'ncircles', 'fg', 'rs', 'bm'
    ], case_sensitive=False),
    default='hamdiv',
    help="Diversity metric to calculate (default: hamdiv)"
)
@click.option(
    "--threshold",
    type=float,
    default=0.7,
    help="Threshold for NCircles metric (default: 0.7)"
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Verbose output"
)
def single(
    smiles_list: tuple[str, ...],
    metric: str,
    threshold: float,
    verbose: bool
):
    """Calculate diversity for a set of SMILES strings provided as arguments.
    
    Examples:
    \b
        hamdiv single "CCO" "CCC" "c1ccccc1"
        hamdiv single "CCO" "CCC" --metric intdiv
        hamdiv single "CCO" "CCC" "c1ccccc1" --metric all
    """
    try:
        # Validate SMILES
        valid_smiles = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                valid_smiles.append(smi)
            elif verbose:
                console.print(f"[yellow]Warning:[/yellow] Invalid SMILES '{smi}' skipped")
        
        if not valid_smiles:
            console.print("[red]Error:[/red] No valid SMILES provided")
            sys.exit(1)
        
        if verbose:
            console.print(f"[green]Processing {len(valid_smiles)} valid molecules[/green]")
        
        # Calculate diversity metrics
        if metric == 'all':
            results = diversity_vector(smiles=valid_smiles)
            # Pretty output for all metrics
            for metric_name, value in results.items():
                if isinstance(value, (int, float)):
                    click.echo(f"{metric_name.replace('_', ' ').title()}: {value:.4f}")
                else:
                    click.echo(f"{metric_name.replace('_', ' ').title()}: {value}")
        else:
            # Calculate specific metric
            if metric == 'ncircles':
                metric_name = f"NCircles-{threshold}"
            else:
                metric_map = {
                    'hamdiv': 'HamDiv',
                    'richness': 'Richness',
                    'intdiv': 'IntDiv',
                    'sumdiv': 'SumDiv',
                    'diam': 'Diam',
                    'sumdiam': 'SumDiam',
                    'bot': 'Bot',
                    'sumbot': 'SumBot',
                    'dpp': 'DPP',
                    'fg': 'FG',
                    'rs': 'RS',
                    'bm': 'BM'
                }
                metric_name = metric_map[metric]
            
            value = diversity_all(smiles=valid_smiles, mode=metric_name)
            click.echo(f"{metric}: {value:.4f}")
        
    except Exception as e:
        console.print(f"[red]Error calculating diversity:[/red] {e}")
        sys.exit(1)


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--format", "input_format",
    type=click.Choice(['smi', 'smiles', 'sdf'], case_sensitive=False),
    help="Input file format (auto-detected if not specified)"
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Verbose output"
)
def validate(
    input_file: Path,
    input_format: Optional[str],
    verbose: bool
):
    """Validate molecules in input file for diversity calculations.
    
    Examples:
    \b
        hamdiv validate molecules.smi
        hamdiv validate structures.sdf --verbose
    """
    if verbose:
        console.print(f"[green]Validating molecules in:[/green] {input_file}")
    
    # Auto-detect format
    if not input_format:
        suffix = input_file.suffix.lower()
        if suffix in ['.smi', '.smiles']:
            input_format = 'smiles'
        elif suffix == '.sdf':
            input_format = 'sdf'
        else:
            input_format = 'smiles'
    
    try:
        smiles_list = read_molecules(input_file, input_format, verbose)
        
        valid_count = 0
        invalid_count = 0
        
        for i, smi in enumerate(smiles_list, 1):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                valid_count += 1
                if verbose:
                    console.print(f"[green]✓[/green] Molecule {i}: {smi}")
            else:
                invalid_count += 1
                if verbose:
                    console.print(f"[red]✗[/red] Molecule {i}: {smi}")
        
        total = valid_count + invalid_count
        console.print(f"\n[yellow]Validation Results:[/yellow]")
        console.print(f"Total molecules: {total}")
        console.print(f"Valid: {valid_count}")
        console.print(f"Invalid: {invalid_count}")
        if total > 0:
            console.print(f"Success rate: {valid_count/total*100:.1f}%")
        
        if invalid_count > 0:
            sys.exit(1)
    
    except Exception as e:
        console.print(f"[red]Error validating molecules:[/red] {e}")
        sys.exit(1)


@cli.command()
def metrics():
    """List available diversity metrics with descriptions."""
    
    metrics_info = {
        'HamDiv': 'Hamiltonian diversity (main novelty - TSP-based)',
        'Richness': 'Number of unique molecules',
        'IntDiv': 'Internal diversity (average pairwise distance)',
        'SumDiv': 'Sum diversity',
        'Diam': 'Diameter (maximum pairwise distance)',
        'SumDiam': 'Sum diameter',
        'Bot': 'Bottleneck (minimum pairwise distance)',
        'SumBot': 'Sum bottleneck',
        'DPP': 'Determinantal point process',
        'NCircles-X': 'N-circles with threshold X (default 0.7)',
        'FG': 'Functional groups diversity',
        'RS': 'Ring systems diversity',
        'BM': 'Bemis-Murcko scaffolds diversity'
    }
    
    table = Table(title="Available Diversity Metrics")
    table.add_column("Metric", style="cyan", no_wrap=True)
    table.add_column("Description", style="white")
    
    for metric, description in metrics_info.items():
        table.add_row(metric, description)
    
    console.print(table)
    
    console.print("\n[yellow]Usage examples:[/yellow]")
    console.print("  hamdiv calculate data.smi --metric hamdiv")
    console.print("  hamdiv calculate data.smi --metric ncircles --threshold 0.8")
    console.print("  hamdiv calculate data.smi --metric all")


@cli.command()
def info():
    """Display system and package information."""
    try:
        from importlib.metadata import version
        hamdiv_version = version('hamdiv')
    except ImportError:
        hamdiv_version = "Unknown"
    
    console.print("[bold]HamDiv Information[/bold]")
    console.print("=" * 30)
    console.print(f"Version: {hamdiv_version}")
    console.print(f"Python: {sys.version.split()[0]}")
    console.print(f"Platform: {sys.platform}")
    console.print("")
    
    # Check dependencies
    deps = {
        'RDKit': 'rdkit',
        'NumPy': 'numpy', 
        'Pandas': 'pandas',
        'NetworkX': 'networkx',
        'Python-TSP': 'python-tsp',
        'tqdm': 'tqdm'
    }
    
    console.print("[yellow]Dependencies:[/yellow]")
    for name, module in deps.items():
        try:
            version_str = version(module)
            console.print(f"  {name}: {version_str}")
        except ImportError:
            console.print(f"  {name}: [red]Not installed[/red]")
    
    console.print("")
    console.print("[yellow]Available metrics:[/yellow]")
    console.print("  13 molecular diversity metrics including Hamiltonian diversity")
    console.print("")
    console.print("[yellow]Features:[/yellow]")
    console.print("  • Hamiltonian diversity (TSP-based)")
    console.print("  • Count-based metrics (Richness, FG, RS, BM)")
    console.print("  • Distance-based metrics (IntDiv, SumDiv, Diam, etc.)")
    console.print("  • ECFP fingerprints with Tanimoto similarity")


def read_molecules(input_file: Path, file_format: str, verbose: bool = False) -> List[str]:
    """Read molecules from input file and return SMILES list."""
    smiles_list = []
    
    if file_format.lower() in ['smi', 'smiles']:
        with open(input_file, 'r') as f:
            lines = f.readlines()
            
        if verbose:
            lines = track(lines, description="Reading SMILES...")
        
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if line and not line.startswith('#'):
                try:
                    # Handle both SMILES and SMILES + name format
                    parts = line.split('\t') if '\t' in line else line.split()
                    smiles = parts[0]
                    smiles_list.append(smiles)
                except Exception as e:
                    if verbose:
                        console.print(f"[yellow]Warning:[/yellow] Error reading line {line_num}: {e}")
    
    elif file_format.lower() == 'sdf':
        try:
            supplier = Chem.SDMolSupplier(str(input_file))
            
            if verbose:
                mol_count = len([mol for mol in Chem.SDMolSupplier(str(input_file)) if mol is not None])
                supplier = Chem.SDMolSupplier(str(input_file))
                mols = list(track(supplier, total=mol_count, description="Reading SDF..."))
            else:
                mols = list(supplier)
            
            # Convert to SMILES
            for mol in mols:
                if mol is not None:
                    try:
                        smi = Chem.MolToSmiles(mol)
                        smiles_list.append(smi)
                    except Exception as e:
                        if verbose:
                            console.print(f"[yellow]Warning:[/yellow] Error converting molecule to SMILES: {e}")
                
        except Exception as e:
            console.print(f"[red]Error reading SDF file:[/red] {e}")
            return []
    
    else:
        console.print(f"[red]Error:[/red] Unsupported format '{file_format}'")
        return []
    
    return smiles_list


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == "__main__":
    main()