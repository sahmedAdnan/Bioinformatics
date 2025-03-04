#!/usr/bin/env python

import random
from common.io import read_fasta, get_substitution_matrix
from common.utils import show_header, loading, validate_sequence, get_int_input, Colors, center_text
from affine.alignment import global_affine_alignment
from linear.alignment import global_linear_alignment

def save_alignment_to_fasta(alignment, file_name):
    try:
        with open(file_name, 'w') as f:
            f.write(">Alignment1\n")
            f.write(alignment[0] + "\n")
            f.write(">Alignment2\n")
            f.write(alignment[1] + "\n")
        print(Colors.OKGREEN + center_text(f"✔ Alignment saved to {file_name}") + Colors.ENDC)
    except Exception as e:
        print(Colors.FAIL + center_text(f"✖ Error saving alignment: {str(e)}") + Colors.ENDC)

def main():
    while True:
        show_header("DNA Sequence Alignment with Gap Cost")
        
        
        # Sequence Input.
        seq1, seq2 = None, None
        print(Colors.OKCYAN + "\n[1] Manual Sequence Input" + Colors.ENDC)
        print(Colors.OKCYAN + "[2] FASTA File Input" + Colors.ENDC)
        print(Colors.OKCYAN + "[3] Exit" + Colors.ENDC)
        choice = input(Colors.WARNING + "\nSelect Sequence Input Method (1/2/3): " + Colors.ENDC).strip()
        if choice == '1':
            while True:
                seq1 = input(Colors.OKBLUE + "\nEnter first sequence (ACGT): " + Colors.ENDC).upper().strip()
                if validate_sequence(seq1):
                    break
            while True:
                seq2 = input(Colors.OKBLUE + "Enter second sequence (ACGT): " + Colors.ENDC).upper().strip()
                if validate_sequence(seq2):
                    break
        elif choice == '2':
            path1 = input(Colors.OKBLUE + "\nEnter first FASTA file path: " + Colors.ENDC).strip()
            path2 = input(Colors.OKBLUE + "Enter second FASTA file path: " + Colors.ENDC).strip()
            seq1 = read_fasta(path1)
            seq2 = read_fasta(path2)
            if not (seq1 and seq2 and validate_sequence(seq1) and validate_sequence(seq2)):
                print(Colors.FAIL + "\nError reading sequences. Please try again." + Colors.ENDC)
                continue
        elif choice == '3':
            print(Colors.OKCYAN + "\nThank you for using our alignment tool. Goodbye!" + Colors.ENDC)
            break
        else:
            print(Colors.FAIL + "\nInvalid choice. Please try again." + Colors.ENDC)
            continue
        
        # Substitution Matrix.
        score_matrix = get_substitution_matrix()
        
        # Gap Cost Method.
        print(Colors.OKCYAN + "\n[1] Global Linear Alignment" + Colors.ENDC)
        print(Colors.OKCYAN + "[2] Global Affine Alignment" + Colors.ENDC)
        while True:
            method_choice = input(Colors.WARNING + "\nSelect Gap Cost Method (1/2): " + Colors.ENDC).strip()
            if method_choice in ['1','2']:
                break
            print(Colors.FAIL + "Invalid choice, please enter 1 or 2." + Colors.ENDC)
        

        if method_choice == '1':
            gap = get_int_input("Enter gap penalty (integer): ", "✖ Please enter a valid integer for gap penalty.")
            loading("Computing Global Linear Alignment...")
            score, alignments = global_linear_alignment(seq1, seq2, score_matrix, gap)
        else:
            go = get_int_input("Enter gap opening penalty (integer): ", "✖ Please enter a valid integer for gap opening.")
            ge = get_int_input("Enter gap extension penalty (integer): ", "✖ Please enter a valid integer for gap extension.")
            loading("Computing Global Affine Alignment...")
            score, alignments = global_affine_alignment(seq1, seq2, score_matrix, go, ge, show_all=True)
        
        # Display Optimal Score.
        print("\n" + Colors.OKCYAN + center_text("Optimal Score: " + str(score)) + Colors.ENDC)
        
        # Ask if user wants to see the alignments.
        show_align = input(Colors.OKBLUE + "\nDisplay the possible alignments? (y/n): " + Colors.ENDC).strip().lower()
        if show_align == 'y':
            print(Colors.BOLD + center_text(f"\nFound {len(alignments)} possible alignment(s)") + Colors.ENDC)
            for idx, (a1, a2) in enumerate(alignments, 1):
                print(Colors.OKGREEN + f"\nAlignment {idx}:" + Colors.ENDC)
                print(Colors.OKBLUE + center_text(a1) + Colors.ENDC)
                print(Colors.OKBLUE + center_text(a2) + Colors.ENDC)
        else:
            print(Colors.OKBLUE + "\nSkipping display of alignments." + Colors.ENDC)
        
        # If multiple alignments exist, ask if user wants to save one randomly.
        if len(alignments) >= 1:
            save_choice = input(Colors.OKBLUE + "\nSave one random alignment to a FASTA file? (y/n): " + Colors.ENDC).strip().lower()
            if save_choice == 'y':
                file_name = input(Colors.OKBLUE + "Enter file name (e.g., alignment.fasta): " + Colors.ENDC).strip()
                chosen_alignment = alignments[random.randint(0, len(alignments)-1)]
                save_alignment_to_fasta(chosen_alignment, file_name)
        

        # Ask to run another alignment.
        again = input(Colors.OKBLUE + "\nPerform another alignment? (y/n): " + Colors.ENDC).strip().lower()
        if again != 'y':
            print(Colors.OKCYAN + "\nThank you for using our alignment tool. Goodbye!" + Colors.ENDC)
            break

if __name__ == "__main__":
    main()
