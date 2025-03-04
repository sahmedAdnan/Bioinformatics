# common/io.py
from Bio import SeqIO
from common.utils import Colors

def read_fasta(file_path):
    try:
        sequences = []
        with open(file_path, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                sequences.append(str(record.seq).upper())
        if not sequences:
            raise ValueError("FASTA file is empty or improperly formatted.")
        print(Colors.OKGREEN + f"✔ Loaded {len(sequences)} sequence(s) from FASTA file." + Colors.ENDC)
        return sequences[0] if len(sequences) == 1 else sequences
    except Exception as e:
        print(Colors.FAIL + f"✖ Error: Could not read FASTA file. {str(e)}" + Colors.ENDC)
        return None

def manual_score_matrix_input():
    print(Colors.HEADER + "\n🔹 Manual Substitution Matrix Entry" + Colors.ENDC)
    try:
        headers = ["A", "C", "G", "T"]
        matrix = {h: {} for h in headers}
        for header in headers:
            while True:
                row = input(f"Enter scores for {header} (space-separated ACGT): ").strip().split()
                if len(row) != 4:
                    print(Colors.FAIL + "✖ Invalid input. Need 4 numbers." + Colors.ENDC)
                    continue
                try:
                    # Convert each input to an integer (supports negatives)
                    scores = [int(score) for score in row]
                except ValueError:
                    print(Colors.FAIL + "✖ Invalid input. Make sure to enter valid integers (negative values allowed)." + Colors.ENDC)
                    continue
                matrix[header] = {h: score for h, score in zip(headers, scores)}
                break
        return matrix
    except Exception as e:
        print(Colors.FAIL + f"✖ Matrix input error: {str(e)}" + Colors.ENDC)
        return None


def load_score_matrix(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        if len(lines) < 5:
            raise ValueError("Invalid matrix format")
        headers = ["A", "C", "G", "T"]
        matrix = {h: {} for h in headers}
        for line in lines[1:5]:
            parts = line.split()
            if parts[0] not in headers or len(parts) != 5:
                raise ValueError("Invalid matrix row")
            row_header = parts[0]
            matrix[row_header] = {h: int(score) for h, score in zip(headers, parts[1:5])}
        return matrix
    except Exception as e:
        print(Colors.FAIL + f"✖ Matrix load error: {str(e)}" + Colors.ENDC)
        return None

def get_substitution_matrix():
    while True:
        print(Colors.WARNING + "\n📌 Choose matrix input method:" + Colors.ENDC)
        print(Colors.OKCYAN + "[1] Manual entry" + Colors.ENDC)
        print(Colors.OKCYAN + "[2] Load from file" + Colors.ENDC)
        choice = input(Colors.OKBLUE + "Select an option (1/2): " + Colors.ENDC).strip()
        if choice == '1':
            matrix = manual_score_matrix_input()
            if matrix:
                return matrix
        elif choice == '2':
            while True:
                path = input(Colors.OKBLUE + "Enter file path: " + Colors.ENDC).strip()
                matrix = load_score_matrix(path)
                if matrix:
                    return matrix
                print(Colors.FAIL + "✖ Invalid file, try again." + Colors.ENDC)
        else:
            print(Colors.FAIL + "✖ Invalid choice." + Colors.ENDC)
