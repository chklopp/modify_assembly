import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
  
  
    def create_sequence_from_clipboard(self, seq_name):
        """Create a novel sequence from clipboard."""
        if self.clipboard is None:
            print("Empty clipboard. Impossible to create a sequence.")
        else :
            self.sequences[seq_name] = SeqRecord(Seq(str(self.clipboard)),id=seq_name,name=seq_name,description=seq_name)
            #print(self.sequences)
            print(f"new sequence {seq_name} added to sequences.")

    
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
            


    def split(self, seq_id, position):
        """Split sequence into two sequences at given position.
        The first part keeps the original name, the second gets -2 suffix."""
        if seq_id in self.sequences:
            seq = self.sequences[seq_id].seq
            # Ensure position is within valid range
            position = max(0, min(position, len(seq)))
            
            # First part (keeps original name)
            first_part = seq[:position]
            # Second part (gets -2 suffix)
            second_part = seq[position:]
            
            # Create new sequence for the second part
            new_seq_id = f"{seq_id}-2"
            self.sequences[new_seq_id] = SeqRecord(Seq(str(second_part)), id=new_seq_id, name=new_seq_id, description=new_seq_id)
            
            # Update original sequence to first part
            self.sequences[seq_id].seq = first_part
            
            print(f"Sequence {seq_id} split at position {position}. Created {new_seq_id}.")
        else:
            print(f"Sequence {seq_id} not found.")


    def filter(self, min_size):
        """Filter sequences by size: keep only sequences with length >= min_size.
        All sequences smaller than min_size are removed from the file."""
        removed = [sid for sid, record in self.sequences.items() if len(record.seq) < min_size]
        for sid in removed:
            del self.sequences[sid]
        print(f"{len(removed)} sequence(s) shorter than {min_size} removed.")
        print(f"{len(self.sequences)} sequence(s) kept.")


    def move(self, seq_id, position):
        """Move a sequence to a new position in the sequence file.
        Position is 1-based (1 = first position).
        The number of sequences is not modified."""
        if seq_id in self.sequences:
            # Convert 1-based position to 0-based index
            pos = max(1, min(position, len(self.sequences)))
            pos -= 1  # Convert to 0-based
            
            # Get all sequence IDs in current order
            all_ids = list(self.sequences.keys())
            
            # Remove the sequence from its current position
            all_ids.remove(seq_id)
            
            # Insert at the desired position (clamp to valid range)
            pos = max(0, min(pos, len(all_ids)))
            all_ids.insert(pos, seq_id)
            
            # Rebuild the sequences dictionary in the new order
            self.sequences = {sid: self.sequences[sid] for sid in all_ids}
            
            print(f"Sequence {seq_id} moved to position {position}.")
        else:
            print(f"Sequence {seq_id} not found.")


    def rename(self, old_id, new_id):
        """Rename a sequence from old_id to new_id.
        The position of the renamed sequence is preserved."""
        if old_id not in self.sequences:
            print(f"Sequence {old_id} not found.")
            return
        
        if new_id in self.sequences:
            print(f"Sequence {new_id} already exists. Renaming aborted.")
            return
        
        # Get current order and position of old_id
        all_ids = list(self.sequences.keys())
        pos = all_ids.index(old_id)
        
        # Get the sequence record
        seq_record = self.sequences[old_id]
        seq_record.id = new_id
        seq_record.name = new_id
        seq_record.description = new_id
        
        # Rebuild dictionary: same order but with new_id at position pos
        new_sequences = {}
        for i, sid in enumerate(all_ids):
            if sid == old_id:
                # This is the position where we insert new_id
                new_sequences[new_id] = seq_record
            else:
                new_sequences[sid] = self.sequences[sid]
        
        self.sequences = new_sequences
        print(f"Sequence renamed from {old_id} to {new_id}.")


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
        elif cmd == "create_sequence_from_clipboard" and len(parts) == 2:
            self.create_sequence_from_clipboard(parts[1])
        elif cmd == "paste" and len(parts) == 3:
            self.paste(parts[1], int(parts[2]))
        elif cmd == "clean" and len(parts) == 2:
            self.clean_sequence(parts[1])
        elif cmd == "split" and len(parts) == 3:
            self.split(parts[1], int(parts[2]))
        elif cmd == "move" and len(parts) == 3:
            self.move(parts[1], int(parts[2]))
        elif cmd == "rename" and len(parts) == 3:
            self.rename(parts[1], parts[2])
        elif cmd == "filter" and len(parts) == 2:
            self.filter(int(parts[1]))
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
