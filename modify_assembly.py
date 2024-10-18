fimport argparse
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
            print(f"{len(records)} séquence(s) chargée(s) depuis {file_path}.")
        except FileNotFoundError:
            print("Fichier non trouvé.")

    def save(self, output_file):
        """Sauvegarde les séquences modifiées dans un fichier FASTA."""
        SeqIO.write(self.sequences.values(), output_file, "fasta")
        print(f"Séquences sauvegardées dans {output_file}.")

    def reverse_complement(self, seq_id):
        """Reverse complément la séquence donnée."""
        if seq_id in self.sequences:
            self.sequences[seq_id].seq = self.sequences[seq_id].seq.reverse_complement()
            print(f"Reverse complément de la séquence {seq_id} exécuté.")
        else:
            print(f"Séquence {seq_id} non trouvée.")

    def remove(self, seq_id):
        """supprime la la séquence donnée."""
        if seq_id in self.sequences:
            del self.sequences[seq_id]
            print(f"suppression de la séquence {seq_id} exécuté.")
        else:
            print(f"Séquence {seq_id} non trouvée.")
            
    def delete(self, seq_id, start, end):
        """Supprime une portion de la séquence."""
        if seq_id in self.sequences:
            seq = self.sequences[seq_id].seq
            new_seq = seq[:start] + seq[end+1:]
            self.sequences[seq_id].seq = new_seq
            print(f"Portion supprimée de {start} à {end} dans {seq_id}.")
        else:
            print(f"Séquence {seq_id} non trouvée.")

    def cut(self, seq_id, start, end):
        """Coupe une portion de la séquence et la stocke dans le presse-papiers."""
        if seq_id in self.sequences:
            seq = self.sequences[seq_id].seq
            self.clipboard = seq[start:end+1]
            new_seq = seq[:start] + seq[end+1:]
            self.sequences[seq_id].seq = new_seq
            print(f"Portion coupée de {start} à {end} dans {seq_id} et stockée.")
        else:
            print(f"Séquence {seq_id} non trouvée.")

    def paste(self, seq_id, position):
        """Colle la séquence du presse-papiers à la position spécifiée dans la séquence cible."""
        if self.clipboard is None:
            print("Le presse-papiers est vide. Impossible de coller.")
            return

        if seq_id in self.sequences:
            seq = self.sequences[seq_id].seq
            new_seq = seq[:position] + self.clipboard + seq[position:]
            self.sequences[seq_id].seq = new_seq
            print(f"Séquence collée à la position {position} dans {seq_id}.")
        else:
            print(f"Séquence {seq_id} non trouvée.")

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
            print(f"Les segments composés uniquement de 'N' ont été ajustés à 100 caractères pour {seq_id}.")
        else:
            print(f"Séquence {seq_id} non trouvée.")
            
    def execute(self, command, output):
        """Analyse et exécute une commande donnée par l'utilisateur."""
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
        elif cmd == "remove" and len(parts) == 2:
            self.remove(parts[1])
        elif cmd == "paste" and len(parts) == 3:
            self.paste(parts[1], int(parts[2]))
        elif cmd == "clean" and len(parts) == 2:
            self.clean_sequence(parts[1])
        else:
            print("Commande non reconnue ou paramètres incorrects.")

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

