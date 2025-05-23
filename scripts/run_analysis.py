#!/usr/bin/env python3
"""
Analysis Runner and Navigation Script
====================================
Quick access script for running analyses from the organized workspace structure.
"""

import os
import sys
from pathlib import Path

def show_workspace_structure():
    """Display the organized workspace structure"""
    print("=" * 80)
    print("📁 ORGANIZED WORKSPACE STRUCTURE")
    print("=" * 80)
    
    structure = {
        "01_original_data": "Original input files (VCF, CSV data)",
        "02_analysis_scripts": "Python scripts and Jupyter notebooks", 
        "03_analysis_results": "Analysis outputs and data tables",
        "04_final_reports": "Documentation and summary reports",
        "05_visualizations": "Generated plots and charts"
    }
    
    for folder, description in structure.items():
        if os.path.exists(f"../{folder}"):
            print(f"📂 {folder}/")
            print(f"   {description}")
            
            # Count files in each folder
            try:
                file_count = len(list(Path(f"../{folder}").rglob("*")))
                print(f"   ({file_count} items)")
            except:
                print("   (access restricted)")
            print()

def show_analysis_options():
    """Show available analysis scripts"""
    print("=" * 80) 
    print("🔬 AVAILABLE ANALYSIS SCRIPTS")
    print("=" * 80)
    
    scripts = {
        "european_novel_variants_discovery.py": "European population-specific novel variants",
        "novel_variants_discovery.py": "General novel variant discovery", 
        "high_impact_variants_analysis.py": "High impact variant identification",
        "chromosome_wise_analysis.py": "Chromosome-by-chromosome analysis",
        "full_analysis_direct.py": "Direct full dataset analysis",
        "complete_analysis.py": "Complete comprehensive analysis",
        "quick_analysis.py": "Quick subset analysis for testing",
        "column_explorer.py": "Data structure exploration"
    }
    
    for i, (script, description) in enumerate(scripts.items(), 1):
        if os.path.exists(script):
            print(f"{i}. {script}")
            print(f"   → {description}")
            print()

def check_data_access():
    """Check if data files are accessible"""
    print("=" * 80)
    print("📊 DATA FILE ACCESS CHECK")
    print("=" * 80)
    
    # File paths - update these to match your data location
    files_to_analyze = [
        "data/your_opencravat_data.csv",  # Replace with your OpenCRAVAT annotated file
        # Add more files here if needed
    ]
    
    for file_path in files_to_analyze:
        if os.path.exists(file_path):
            file_size = os.path.getsize(file_path) / (1024*1024*1024)  # GB
            print(f"✅ {os.path.basename(file_path)} ({file_size:.1f} GB)")
        else:
            print(f"❌ {os.path.basename(file_path)} - NOT FOUND")
    print()

def show_key_results():
    """Display key analysis results summary"""
    print("=" * 80)
    print("📈 KEY ANALYSIS RESULTS SUMMARY")
    print("=" * 80)
    
    print("🧬 European Novel Variants Analysis:")
    print("   • 20,128 European novel variants identified")
    print("   • 162 population-specific variants")
    print("   • 1,365 ultra-rare European variants")
    print("   • Top variant: TYK2 (score: 78.2)")
    print()
    
    print("⚡ High Impact Variants Analysis:")
    print("   • 31,105 high impact variants (2.25% of total)")
    print("   • 104 pathogenic variants from ClinVar")
    print("   • 1,982 ultra-high CADD score variants")
    print("   • 1,077 loss-of-function variants")
    print()
    
    print("🧮 Dataset Overview:")
    print("   • Total variants: 1,382,711")
    print("   • Data size: 1.5GB")
    print("   • Annotation columns: 121")
    print("   • Processing: Chunk-based for efficiency")
    print()

def run_analysis(script_name):
    """Run a specific analysis script"""
    if os.path.exists(script_name):
        print(f"🚀 Running {script_name}...")
        print("=" * 80)
        os.system(f"python {script_name}")
    else:
        print(f"❌ Script {script_name} not found!")

def main():
    """Main navigation function"""
    print("🧬 OpenCRAVAT Analysis Workspace Navigator")
    print("Organized structure for efficient variant analysis")
    print()
    
    while True:
        print("\nChoose an option:")
        print("1. Show workspace structure")
        print("2. Show available analysis scripts") 
        print("3. Check data file access")
        print("4. Show key results summary")
        print("5. Run specific analysis")
        print("6. Exit")
        
        choice = input("\nEnter choice (1-6): ").strip()
        
        if choice == "1":
            show_workspace_structure()
        elif choice == "2":
            show_analysis_options()
        elif choice == "3":
            check_data_access()
        elif choice == "4":
            show_key_results()
        elif choice == "5":
            show_analysis_options()
            script_choice = input("\nEnter script name to run: ").strip()
            run_analysis(script_choice)
        elif choice == "6":
            print("👋 Goodbye!")
            break
        else:
            print("❌ Invalid choice. Please enter 1-6.")

if __name__ == "__main__":
    main() 