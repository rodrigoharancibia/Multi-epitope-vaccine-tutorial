import os
import glob
import itertools
import string
import math
import random
from datetime import datetime

# --- CONFIGURATION ---
CONFIDENCE_LEVEL_Z = 1.96   # 95% Confidence


def get_clean_txt_files():
    """Returns .txt files excluding system/output files."""
    files = glob.glob("*.txt")
    ignored = ["map.txt", "requirements.txt", "log.txt"]
    return [f for f in files if f not in ignored]


def select_input_file():
    """Auto-detects single file or prompts user."""
    files = get_clean_txt_files()

    if not files:
        print("\n[ERROR] No valid .txt input file found.")
        return None

    if len(files) == 1:
        print(f"\nInput file detected: '{files[0]}'")
        return files[0]

    print("\nFiles found:")
    for i, f in enumerate(files):
        print(f"  [{i+1}] {f}")

    while True:
        try:
            choice = input("Select file number: ")
            if choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(files):
                    return files[idx]
        except ValueError:
            pass


def parse_epitopes(file_path):
    """Reads epitopes and maps them to labels."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            epitopes = [line.strip() for line in f if line.strip()]

        if not epitopes:
            return None, None

        labels = list(string.ascii_uppercase)
        if len(epitopes) > len(labels):
            labels = [str(i+1) for i in range(len(epitopes))]
        else:
            labels = labels[:len(epitopes)]

        return epitopes, {labels[i]: epitopes[i] for i in range(len(epitopes))}
    except Exception as e:
        print(f"File Read Error: {e}")
        return None, None


def get_linker_selection():
    """Displays ASCII menu and returns selected linker."""
    menu = """

    
┌──────────────────────────┐
│   Select Linker          │
├──────────────────────────┤
│ [1] GPGPG                │
│ [2] KKK                  │
│ [3] GSGSG                │
│ [4] EAAAK                │
│ [5] AAY                  │
│                          │
│ [0] Manual Input         │
└──────────────────────────┘
    """
    print(menu)

    mapping = {
        "1": "GPGPG", "2": "KKK", "3": "GSGSG",
        "4": "EAAAK", "5": "AAY"
    }

    while True:
        opt = input(">> Option: ").strip()
        if opt == "0":
            return input("Enter custom linker sequence: ").strip().upper()
        if opt in mapping:
            return mapping[opt]


def get_sampling_strategy():
    """Displays sampling strategy menu."""
    menu = """

┌─────────────────────────────────────────────────┐
│   SAMPLING STRATEGY                             │
├─────────────────────────────────────────────────┤
│ [1] Latin Square Only                           │
│     • n models (minimal)                        │
│     • 100% positional coverage                  │
│     • ~n unique pairs covered                   │
│                                                 │
│ [2] Latin Square + All Adjacent Pairs           │
│     • ~2n to 3n models (recommended)            │
│     • 100% positional coverage                  │
│     • 100% pairwise interactions                │
│                                                 │
│ [3] Random Sampling (Monte Carlo)               │
│     • Custom sample size                        │
│     • Statistical representation                │
│                                                 │
│ [4] Generate ALL permutations                   │
│     • n! models (only for small n)              │
│                                                 │
└─────────────────────────────────────────────────┘
    """
    print(menu)

    while True:
        opt = input(">> Strategy: ").strip()
        if opt in ["1", "2", "3", "4"]:
            return opt


def generate_latin_square(epitope_keys):
    """Generates Latin Square using circular rotation."""
    n = len(epitope_keys)
    models = []
    
    for k in range(n):
        model = []
        for j in range(n):
            index = (k + j) % n
            model.append(epitope_keys[index])
        models.append(tuple(model))
    
    return models


def extract_pairs_from_models(models):
    """Extracts all adjacent pairs from list of models."""
    pairs = set()
    for model in models:
        for i in range(len(model) - 1):
            pair = (model[i], model[i+1])
            pairs.add(pair)
    return pairs


def generate_models_for_missing_pairs(existing_models, epitope_keys, missing_pairs):
    """
    Greedy algorithm to generate additional models that cover missing pairs.
    Returns list of new models.
    """
    additional_models = []
    remaining_pairs = missing_pairs.copy()
    
    # Generate all possible permutations
    all_permutations = list(itertools.permutations(epitope_keys))
    
    # Filter out already used permutations
    existing_set = set(existing_models)
    available_perms = [p for p in all_permutations if p not in existing_set]
    
    iteration = 0
    max_iterations = len(available_perms)
    
    while remaining_pairs and iteration < max_iterations:
        best_model = None
        best_coverage = 0
        
        # Find permutation that covers most remaining pairs
        for perm in available_perms:
            pairs_in_model = set()
            for i in range(len(perm) - 1):
                pair = (perm[i], perm[i+1])
                if pair in remaining_pairs:
                    pairs_in_model.add(pair)
            
            if len(pairs_in_model) > best_coverage:
                best_coverage = len(pairs_in_model)
                best_model = perm
        
        if best_model and best_coverage > 0:
            additional_models.append(best_model)
            available_perms.remove(best_model)
            
            # Remove covered pairs
            for i in range(len(best_model) - 1):
                pair = (best_model[i], best_model[i+1])
                remaining_pairs.discard(pair)
        else:
            break
        
        iteration += 1
    
    return additional_models


def calculate_pair_coverage(models, n):
    """Calculates statistics about pair coverage."""
    covered_pairs = extract_pairs_from_models(models)
    total_possible_pairs = n * (n - 1)  # Directed pairs
    
    return {
        'covered': len(covered_pairs),
        'total': total_possible_pairs,
        'percentage': (len(covered_pairs) / total_possible_pairs * 100) if total_possible_pairs > 0 else 0,
        'pairs': covered_pairs
    }


def calculate_margin_of_error(n, N):
    """Calculates Margin of Error (E)."""
    if n >= N:
        return 0.0
    p = 0.5
    standard_error = math.sqrt((p * (1 - p)) / n)
    fpc = math.sqrt((N - n) / (N - 1)) if N > 1 else 1
    return CONFIDENCE_LEVEL_Z * standard_error * fpc * 100


def display_stats_dashboard(N, n, error, strategy, pair_stats=None):
    """Renders the statistical analysis ASCII box."""
    N_fmt = f"{N:,}".replace(",", ".")
    n_fmt = f"{n:,}".replace(",", ".")
    err_fmt = f"{error:.4f}%"
    pct_fmt = f"{(n/N*100):.4f}%" if N > 0 else "N/A"

    print(f"\n╔══════════════════════════════════════════════════════╗")
    print(f"║   STATISTICAL SAMPLING ANALYSIS                      ║")
    print(f"╠══════════════════════════════════════════════════════╣")
    print(f"║ Strategy             : {strategy:<30}║")
    print(f"║ Total Population (N) : {N_fmt:<30}║")
    print(f"║ Generated Sample (n) : {n_fmt:<30}║")
    print(f"║ % Represented        : {pct_fmt:<30}║")
    print(f"║ Confidence Level (Z) : 95%                           ║")
    print(f"║ Margin of Error (E)  : {err_fmt:<30}║")
    
    if pair_stats:
        pair_pct = f"{pair_stats['percentage']:.2f}%"
        pair_cov = f"{pair_stats['covered']}/{pair_stats['total']}"
        print(f"╠══════════════════════════════════════════════════════╣")
        print(f"║ PAIRWISE INTERACTION COVERAGE                        ║")
        print(f"╠══════════════════════════════════════════════════════╣")
        print(f"║ Adjacent Pairs       : {pair_cov:<30}║")
        print(f"║ Coverage             : {pair_pct:<30}║")
    
    print(f"╚══════════════════════════════════════════════════════╝\n")


def main():
    print("\n╔════════════════════════╗")
    print("║ EPITOPE RECOMBINATOR   ║")
    print("╚════════════════════════╝")

    # 1. Input
    fpath = select_input_file()
    if not fpath:
        return

    ep_list, ep_map = parse_epitopes(fpath)
    if not ep_list:
        return

    # 2. Config
    linker = get_linker_selection()
    strategy = get_sampling_strategy()
    
    k = len(ep_list)
    ep_keys = list(ep_map.keys())
    total_pop = math.perm(k, k)

    # 3. Generate models based on strategy
    permutations = []
    strategy_name = ""
    
    if strategy == "1":
        # Latin Square Only
        strategy_name = "Latin Square (Positional Coverage)"
        permutations = generate_latin_square(ep_keys)
        
    elif strategy == "2":
        # Latin Square + All Adjacent Pairs
        strategy_name = "Latin Square + Pairwise Coverage"
        
        # Step 1: Generate Latin Square
        ls_models = generate_latin_square(ep_keys)
        
        # Step 2: Identify missing pairs
        covered_pairs = extract_pairs_from_models(ls_models)
        all_possible_pairs = set()
        for e1 in ep_keys:
            for e2 in ep_keys:
                if e1 != e2:
                    all_possible_pairs.add((e1, e2))
        
        missing_pairs = all_possible_pairs - covered_pairs
        
        # Step 3: Generate additional models
        additional_models = generate_models_for_missing_pairs(ls_models, ep_keys, missing_pairs)
        
        permutations = ls_models + additional_models
        
    elif strategy == "3":
        # Monte Carlo Random Sampling
        strategy_name = "Random Sampling (Monte Carlo)"
        
        print(f"\nTotal possible permutations: {total_pop:,}")
        while True:
            try:
                sample_size = int(input("Enter sample size (or 0 for auto): "))
                if sample_size == 0:
                    # Auto-calculate using sqrt(N) rule
                    sample_size = min(int(math.sqrt(total_pop)) + k, total_pop)
                    print(f"Auto-calculated sample size: {sample_size}")
                    break
                elif 0 < sample_size <= total_pop:
                    break
                else:
                    print(f"Please enter a value between 1 and {total_pop}")
            except ValueError:
                print("Invalid input. Please enter a number.")
        
        if sample_size >= total_pop:
            permutations = list(itertools.permutations(ep_keys))
        else:
            seen = set()
            while len(permutations) < sample_size:
                shuffled = random.sample(ep_keys, k)
                perm_tuple = tuple(shuffled)
                if perm_tuple not in seen:
                    seen.add(perm_tuple)
                    permutations.append(perm_tuple)
    
    elif strategy == "4":
        # All permutations
        strategy_name = "Complete Enumeration (All Permutations)"
        
        if total_pop > 50000:
            confirm = input(f"\n⚠️  WARNING: This will generate {total_pop:,} files. Continue? (yes/no): ")
            if confirm.lower() != "yes":
                print("Operation cancelled.")
                return
        
        permutations = list(itertools.permutations(ep_keys))

    # 4. Calculate statistics
    sample_size = len(permutations)
    margin_error = calculate_margin_of_error(sample_size, total_pop)
    pair_stats = calculate_pair_coverage(permutations, k)
    
    display_stats_dashboard(total_pop, sample_size, margin_error, strategy_name, pair_stats)

    # 5. Generate output directory
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_dir = f"epitope_combinations_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)

    # 6. Audit
    used_keys = set(key for perm in permutations for key in perm)
    missing = set(ep_keys) - used_keys
    
    if missing:
        print(f"[Audit] WARNING: Missing coverage for epitopes: {missing}")

    # 7. Write FASTA files
    map_log = []

    for i, perm in enumerate(permutations):
        seq_parts = []
        desc_parts = []

        # Add epitopes with linkers
        for ep_key in perm[:-1]:
            ep_seq = ep_map[ep_key]
            seq_parts.extend([ep_seq, linker])
            desc_parts.extend([f"Ep_{ep_key}", linker])

        # Add last epitope + His-tag
        last_key = perm[-1]
        seq_parts.extend([ep_map[last_key], "HHHHHH"])
        desc_parts.extend([f"Ep_{last_key}", "6xHis"])

        full_seq = "".join(seq_parts)
        full_desc = "-".join(desc_parts)

        filename = f"Construct_{i+1:04d}.fasta"
        header = f">Construct_{i+1} | {full_desc}"

        with open(os.path.join(output_dir, filename), "w") as f:
            f.write(f"{header}\n{full_seq}\n")

        map_log.append(f"{filename}\t{full_desc}")

    # 8. Write map file with detailed statistics
    with open(os.path.join(output_dir, "map.txt"), "w") as f:
        f.write("=" * 70 + "\n")
        f.write("CHIMERIC PROTEIN GENERATION REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Generated: {timestamp}\n")
        f.write(f"Strategy: {strategy_name}\n\n")
        
        f.write("EPITOPE MAPPING:\n")
        f.write("-" * 70 + "\n")
        for key in sorted(ep_map.keys()):
            f.write(f"  Ep_{key}: {ep_map[key]}\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("STATISTICAL ANALYSIS:\n")
        f.write("=" * 70 + "\n")
        f.write(f"  Total Permutations (N): {total_pop:,}\n")
        f.write(f"  Generated Models (n):   {sample_size:,}\n")
        f.write(f"  Percentage Sampled:     {(sample_size/total_pop*100):.4f}%\n")
        f.write(f"  Confidence Level:       95%\n")
        f.write(f"  Margin of Error:        {margin_error:.4f}%\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("PAIRWISE COVERAGE ANALYSIS:\n")
        f.write("=" * 70 + "\n")
        f.write(f"  Adjacent Pairs Covered: {pair_stats['covered']}/{pair_stats['total']}\n")
        f.write(f"  Coverage Percentage:    {pair_stats['percentage']:.2f}%\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("CONSTRUCT DETAILS:\n")
        f.write("=" * 70 + "\n")
        f.write("Filename\t\tStructure\n")
        f.write("-" * 70 + "\n")
        f.write("\n".join(map_log))
        
        # Add coverage matrix
        f.write("\n\n" + "=" * 70 + "\n")
        f.write("POSITIONAL COVERAGE MATRIX:\n")
        f.write("=" * 70 + "\n")
        f.write("Shows which epitope appears in which position\n\n")
        
        # Build matrix
        position_coverage = {ep: set() for ep in ep_keys}
        for perm in permutations:
            for pos_idx, ep_key in enumerate(perm):
                position_coverage[ep_key].add(pos_idx + 1)
        
        f.write(f"{'Epitope':<10}")
        for pos in range(1, k+1):
            f.write(f"Pos{pos:<4}")
        f.write("\n" + "-" * 70 + "\n")
        
        for ep_key in sorted(ep_keys):
            f.write(f"Ep_{ep_key:<7}")
            for pos in range(1, k+1):
                mark = "✓" if pos in position_coverage[ep_key] else "✗"
                f.write(f"{mark:<5}")
            f.write("\n")

    # 9. Final output
    print(f"✓ SUCCESS! {len(permutations)} constructs generated.")
    print(f"✓ LOCATION: {os.path.abspath(output_dir)}")
    print(f"✓ Pairwise coverage: {pair_stats['percentage']:.2f}%")


if __name__ == "__main__":
    main()


