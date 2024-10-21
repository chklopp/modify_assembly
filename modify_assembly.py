import argparse
from Bio import SeqIO
from Bio.Seq import Seq

class SequenceManipulator:
    def __init__(self):
        self.sequences = {}  # Dictionnaire pour stocker les séquences chargées
        self.clipboard = None  # Espace temporaire pour les opérations "couper/coller"

    def load(self, file_path):
        """Charge les séquences d'un fichier FASTA."""
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
            for record in records:
                self.sequences[record.id] = record
            print(f"{len(records)} sequences loaded from {file_path}.")
        except FileNotFoundError:
            print("File not found.")

    def save(self, output_file):
        """Saving sequences from a modified FASTA."""
        SeqIO.write(self.sequences.values(), output_file, "fasta")
        print(f"Sequences saved in {output_file}.")

    def reverse_complement(self, seq_id):
        """Reverse complement the given sequence."""
        if seq_id in self.sequences:
            self.sequences[seq_id].seq = self.sequences[seq_id].seq.reverse_complement()
            print(f"Reverse complément de la séquence {seq_id} exécuté.")
        else:
            print(f"Sequence {seq_id} not found.")

    def remove(self, seq_id):
        """remove sequence from dictionnary."""
        if seq_id in self.sequences:
            del self.sequences[seq_id]
            print(f"suppression de la séquence {seq_id} exécuté.")
        else:
            print(f"Sequence {seq_id} not found.")
            
    def delete(self, seq_id, start, end):
        """Remove sequence protion."""
        if seq_id in self.sequences:
            seq = self.sequences[seq_id].seq
            new_seq = seq[:start] + seq[end+1:]
            self.sequences[seq_id].seq = new_seq
            print(f"Portion supprimée de {start} à {end} dans {seq_id}.")
        else:
            print(f"Sequence {seq_id} not found.")

    def cut(self, seq_id, start, end):
        """Cut sequence and stores the cut part in clipboard."""
        if seq_id in self.sequences:
            seq = self.sequences[seq_id].seq
            self.clipboard = seq[start:end+1]
            new_seq = seq[:start] + seq[end+1:]
            self.sequences[seq_id].seq = new_seq
            print(f"Portion cut {start} à {end} in {seq_id} and stored in clipboard.")
        else:
            print(f"Sequence {seq_id} not found.")
    
    def reverse_complement_clipboard(self):
        """Reverse complement clipboard."""
        if self.clipboard is None:
            print("Empty clipboard. Impossible to reverse it.")
        else :
            self.clipboard = self.clipboard.reverse_complement()
            print("Clipboad reversed.")

    def paste(self, seq_id, position):
        """Paste clipboard at specific location in a given sequence."""
        if self.clipboard is None:
            print("Empty clipboard. Impossible to paste it.")
            return

        if seq_id in self.sequences:
            seq = self.sequences[seq_id].seq
            new_seq = seq[:position] + self.clipboard + seq[position:]
            self.sequences[seq_id].seq = new_seq
            print(f"Sequence pasted at position {position} in {seq_id}.")
        else:
            print(f"Sequence {seq_id} not found.")

    def clean_sequence(self, seq_id):
        """Recherche les segments composés uniquement de 'N' et ajuste leur longueur à 100 'N'."""
        if seq_id in self.sequences:
            seq = str(self.sequences[seq_id].seq)
            new_seq = ""
            i = 0
            while i < len(seq):
                if seq[i] == "N":
                    # Trouver la fin de la séquence de 'N'
                    start = i
                    while i < len(seq) and seq[i] == "N":
                        i += 1
                    # Remplacer la séquence par 100 'N'
                    new_seq += "N" * 100
                else:
                    new_seq += seq[i]
                    i += 1
            self.sequences[seq_id].seq = Seq(new_seq)
            print(f"Segments made of 'N's have been adjusted to 100 characters for sequence {seq_id}.")
        else:
            print(f"Sequence {seq_id} not found.")
            
    def execute(self, command, output):
        """Analyse et execute user commands."""
        parts = command.split()
        if not parts:
            print("Commande vide.")
            return

        cmd = parts[0]
        if cmd == "load" and len(parts) == 2:
            self.load(parts[1])
        elif cmd == "save" and len(parts) == 1:
            self.save(output)
        elif cmd == "reverse_complement" and len(parts) == 2:
            self.reverse_complement(parts[1])
        elif cmd == "delete" and len(parts) == 4:
            self.delete(parts[1], int(parts[2]), int(parts[3]))
        elif cmd == "cut" and len(parts) == 4:
            self.cut(parts[1], int(parts[2]), int(parts[3]))
        elif cmd == "reverse_complement_clipboard" and len(parts) == 1:
            self.reverse_complement_clipboard()
        elif cmd == "remove" and len(parts) == 2:
            self.remove(parts[1])
        elif cmd == "paste" and len(parts) == 3:
            self.paste(parts[1], int(parts[2]))
        elif cmd == "clean" and len(parts) == 2:
            self.clean_sequence(parts[1])
        else:
            print("Non recognized command or wrong number of parameters.")

def main():
    parser = argparse.ArgumentParser(description="Manipulation de fichiers FASTA via commandes.")
    parser.add_argument("--input", type=str, required=True, help="Chemin du fichier FASTA à charger.")
    parser.add_argument("--commands", type=str, required=True, help="Chemin du fichier contenant les commandes.")
    parser.add_argument("--output", type=str, required=True, help="Chemin du fichier de sortie.")

    args = parser.parse_args()

    # Créer une instance du manipulateur de séquences
    manipulator = SequenceManipulator()

    # Charger le fichier FASTA initial
    manipulator.load(args.input)

    # Lire les commandes du fichier de commandes
    with open(args.commands, "r") as cmd_file:
        for line in cmd_file:
            command = line.strip()
            if command:
                manipulator.execute(command, args.output)

    # Sauvegarder les résultats dans le fichier de sortie
    manipulator.save(args.output)

if __name__ == "__main__":
    main()

