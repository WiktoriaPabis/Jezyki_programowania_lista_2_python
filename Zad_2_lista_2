class DNASequence:
    """
    Klasa reprezentująca sekwencję DNA.

    Attributes:
        identifier: Unikalny identyfikator sekwencji.
        data: Sekwencja DNA.
        length1: Długość sekwencji DNA.
        valid_Chars: Zbiór prawidłowych znaków dla sekwencji DNA ('A', 'T', 'C', 'G').
    
    Args:
        number: Liczba porządkowa sekwencji.
        DNA: Sekwencja znaków nici kodujących DNA (złożona z liter 'A', 'T', 'C', 'G').

    Raises:
        ValueError: W przypadku, gdy sekwencja DNA zawiera inne znaki niż 'A', 'T', 'C', 'G'.
    """
    def __init__(self, number, DNA):
        self.identifier = number
        self.valid_Chars = {'A', 'T', 'C', 'G'}
        
        for elem in DNA:
            if elem not in self.valid_Chars:
                raise ValueError("Nieprawidłowe wartości dla sekwencji DNA")
        self.data = DNA
        self.length1 = len(DNA)

    def DNAtext(self):
        """
        Metoda słuząca do wyświetlania wartości identifier oraz nici DNA.

        Returns:
            Tekstowa reprezentacja sekwencji DNA.
        """
        Text = f"> {self.identifier}\n {self.data}"
        return Text

    def mutate(self, position, value):
        """
        Metoda zamienia zasadę na zadanej pozycji sekwencji DNA przez zastąpienie znaku na określonej pozycji.

        Args:
            position: Pozycja do zmiany w sekwencji DNA.
            value: Nowa wartość, która zostanie wpisana w miejsce konkretnej pozycji.

        Raises:
            ValueError: W przypadku, gdy podany jest nieprawidłowy znak oraz, gdy liczba wykracza poza zakres nici DNA.
        """
        dataTab = list(self.data)
        
        if value not in self.valid_Chars:
            raise ValueError("Nieprawidłowa wartość dla sekwencji DNA")
        if len(dataTab) < position:
            raise ValueError("Liczba wykraczająca poza zakres nici DNA")
        dataTab[position - 1] = value
        # join() łączy wszystkie elementy listy w ciąg znaków
        self.data = ''.join(dataTab) 

    def findMotif(self, motif):
        """
        Metoda znajduje wszystkie pozycje, na których występuje zadany motyw w sekwencji DNA.

        Args:
            motif: Motyw, który moa zostać znaleziony w sekwencji.

        Returns:
            Lista pozycji, gdzie motyw występuje w sekwencji DNA.

        Raises:
            ValueError: W przypadku, gdy podana wartość jest nieprawidłowa dla sekwencji DNA.
        """
        if motif not in self.valid_Chars:
            raise ValueError("Nieprawidłowe wartości dla sekwencji DNA")
    
        listOfNucleotides = []
        for index, elem in enumerate(self.data):
            if elem == motif:
                listOfNucleotides.append(index + 1)
    
        return listOfNucleotides
    
    def Complement(self):
        """
        Metoda zmienia sekwencję nici kodującej DNA w kierunku od 5' do 3' na jej komplementarną sekwencję nici matrycowej DNA w kierunku 
        od 3' do 5'.

        Returns:
            Komplementarna matrycowa sekwencja DNA.

        Raises:
            ValueError: W przypadku, gdy podana zasada nie istnieje.
        """
        DNA_mat = ""
        for nucleotide in self.data:
            if nucleotide == 'A':
                DNA_mat = DNA_mat + 'T'
            elif nucleotide == 'T':
                DNA_mat = DNA_mat + 'A'
            elif nucleotide == 'C':
                DNA_mat = DNA_mat + 'G'
            elif nucleotide == 'G':
                DNA_mat = DNA_mat + 'C'
            else:
                raise ValueError("Podano nieprawidłowe wartości")

        DNA_mat = DNA_mat[::-1] #odwrócenie nici

        return DNA_mat

    def transcribe(self, DNA_matryc):
        """
        Metoda przyjmuje sekwencję nici matrycowej DNA w kierunku od 3' do 5' i zwraca jej komplementarną sekwencję RNA w kierunku od 5' do 3'..

        Args:
            DNA_matryc: Ciąg znaków rezprezentujący sekwencję nici metrycowej DNA.

        Returns:
            Obiekt reprezentujący sekwencję nici RNA.

        Raises:
            ValueError: W przypadku, gdy podana zasada nie istnieje.
        """
        RNA = ""
        DNA_matryc_1 = DNA_matryc[::-1]  #odwrócna nić

        for zasada in DNA_matryc_1:
            if zasada == 'A':
                RNA = RNA + 'U'
            elif zasada == 'T':
                RNA = RNA + 'A'
            elif zasada == 'C':
                RNA = RNA + 'G'
            elif zasada == 'G':
                RNA = RNA + 'C'
            else:
                raise ValueError("Podano nieprawidłowe wartości")

        return RNASequence(self.identifier, RNA)
    
class RNASequence:
    """
    Klasa reprezentująca sekwencję RNA.

    Attributes:
        identifier: Unikalny identyfikator sekwencji.
        data: Sekwencja RNA.
        length1: Długość sekwencji RNA.
        valid_Chars: Zbiór prawidłowych znaków dla sekwencji RNA ('A', 'U', 'C', 'G').
        codonsDict: Słownik kodonów RNA i odpowiadających im aminokwasów.

    Args:
        number: Liczba porządkowa sekwencji.
        RNA: Sekwencja RNA.

    Raises:
        ValueError: W przypadku, gdy podane znaki są niedozwolone w nici RNA.
    """
    def __init__(self, number, RNA):
        self.identifier = number
        self.valid_Chars = {'A', 'U', 'C', 'G'}
        for elem in RNA:
            if elem not in self.valid_Chars:
                raise ValueError("Nieprawidłowe wartości dla sekwencji RNA")
        self.data = RNA
        self.length1 = len(RNA)
        self.codonsDict = self.codons()

    def codons(self):
        """
        Metoda zwracająca słownik kodonów RNA na aminokwasy.

        Returns:
            Słownik zawierający kodony RNA i odpowiadające im aminokwasy.

        Słownik wygenerowany za pomocą modelu językowego.
        """
        codonsDict = {
            # Alanina
            "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
            # Cysteina
            "UGU": "Cys", "UGC": "Cys",
            # Kwas asparaginowy
            "GAU": "Asp", "GAC": "Asp",
            # Kwas glutaminowy
            "GAA": "Glu", "GAG": "Glu",
            # Fenyloalanina
            "UUU": "Phe", "UUC": "Phe",
            # Glicyna
            "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
            # Histydyna
            "CAU": "His", "CAC": "His",
            # Izoleucyna
            "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",
            # Lizyna
            "AAA": "Lys", "AAG": "Lys",
            # Leucyna
            "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
            # Metionina lub kodon rozpoczynający transkrypcję
            "AUG": "Met",
            # Asparagina
            "AAU": "Asn", "AAC": "Asn",
            # Prolina
            "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
            # Glutamina
            "CAA": "Gln", "CAG": "Gln",
            # Arginina
            "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
            # Seryna
            "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser",
            # Treonina
            "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
            # Walina
            "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
            # Tryptofan
            "UGG": "Trp",
            # Tyrozyna
            "UAU": "Tyr", "UAC": "Tyr",
            # Kodony STOP
            "UAA": "Ter", "UAG": "Ter", "UGA": "Ter"
        }
    
        return codonsDict

    def RNAtext(self):
        """
        Metoda zwraca tekstową reprezentację sekwencji RNA.

        Returns:
            Tekstowa reprezentacja sekwencji RNA.
        """
        Text1 = f"> {self.identifier}\n {self.data}"

        return Text1
    
    def mutate(self, position, value):
        """
        Metoda zamienia zasadę na zadanej pozycji sekwencji RNA przez zastąpienie znaku na określonej pozycji.

        Args:
            position: Pozycja do zmiany w sekwencji RNA.
            value: Nowa wartość, która zostanie wpisana w miejsce konkretnej pozycji.

        Raises:
            ValueError: W przypadku, gdy podany jest nieprawidłowy znak oraz, gdy liczba wykracza poza zakres nici RNA.
        """
        dataTab = list(self.data)
        if value not in self.valid_Chars:
            raise ValueError("Nieprawidłowa wartość dla sekwencji RNA")
        if len(dataTab) < position:
            raise ValueError("Liczba wykraczająca poza zakres nici RNA")
        dataTab[position - 1] = value
        self.data = ''.join(dataTab)

    def findMotif(self, motif):
        """
        Metoda znajduje wszystkie pozycje, na których występuje zadany motyw w sekwencji RNA, motyw składa się z trzech znaków.

        Args:
            motif: Motyw, który moa zostać znaleziony w sekwencji.

        Returns:
            Lista pozycji, gdzie motyw występuje w sekwencji RNA.

        Raises:
            ValueError: W przypadku, gdy podano nieprawidłową ilość zasad azotowych oraz, gdy podano nieprawidłową wartość dla sekwencji RNA.
        """
        if motif[0] not in self.valid_Chars or motif[1] not in self.valid_Chars or motif[2] not in self.valid_Chars:
            raise ValueError("Nieprawidłowe wartości dla sekwencji RNA")
    
        listOfNucleotides = []
        counter = 0
    
        if len(self.data) % 3 != 0:
            raise ValueError("Nieprawidłowa liczba zasad azotowych, proszę o wprowadzenie ilość zasad azotowych podzielną przez liczbę 3")
    
        for index, elem in enumerate(self.data):
            if counter == 2:
                if elem == motif[2] and self.data[index - 1] == motif[1] and self.data[index - 2] == motif[0]:
                    listOfNucleotides.append((index + 1) // 3) # zabezpieczenie- dzielenie całkowitoliczbowe
                counter = -1
            counter += 1
    
        return listOfNucleotides
    
    def transcribe(self):
        """
        Metoda przekształca sekwencję RNA na sekwencję białkową.

        Returns:
            Obiekt klasy ProteinSequence.

        Raises:
            ValueError: W przypadku, gdy liczba zasad azotowych nie jest podzielma przez 3, lub, gdy nie występuje kodon STOP lub rozpoczynający.
        """
        listAmino = []
        tempCodon = ""
        startTransl = False
        stopTransl = False
        tempListOfCodons = []
        count = 0
    
        if len(self.data) % 3 != 0:
            raise ValueError("Nieprawidłowa liczba zasad azotowych, proszę o wprowadzenie ilość zasad azotowych podzielną przez liczbę 3")
    
        for elem in self.data:
            tempCodon = tempCodon + elem
            if count == 2:
                if startTransl:
                    # jesli translacja rozpoczeta i jest jeden z kodonow STOP
                    if tempCodon == "UAA" or tempCodon == "UAG" or tempCodon == "UGA":
                        stopTransl = True
                        tempListOfCodons.append(tempCodon)
                if tempCodon == "AUG" or startTransl and not stopTransl:
                    tempListOfCodons.append(tempCodon)
                    startTransl = True
            
                tempCodon = ""
                # chce znowu miec count 0 w kolejnej iteracji petli
                count = -1
            count += 1
    
        for elem in tempListOfCodons:
            # elem to klucz, chce dodac value czyli np. Ala
            listAmino.append(self.codonsDict[elem])
    
        if not startTransl:
            raise ValueError("Brak kodonu rozpoczynającego")
        if not stopTransl:
            raise ValueError("Brak kodonu STOP")
    
        return ProteinSequence(self.identifier, listAmino)
    
class ProteinSequence:
    """
    Klasa reprezentująca sekwencję protein.

    Attributes:
        identifier: Unikalny identyfikator sekwencji.
        amino: Lista aminokwasów reprezentujących sekwencję protein.
        length1: Długość sekwencji protein.
        valid_Chars: Zbiór prawidłowych znaków dla sekwencji protein.
        aminoFASTADict: Słownik aminokwasów potrzebny do stworzenia formatu FASTA.
    
    Args:
        number: Liczba porządkowa sekwencji.
        protein: Lista znaków protein.

    Raises:
        ValueError: W przypadku, gdy podane znaki sa niedozwolone dla protein.
    """
    def __init__(self, number, protein):
        self.identifier = number
        self.valid_Chars = {"Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
                            "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "Ter"}
        for elem in protein:
            if elem not in self.valid_Chars:
                raise ValueError("Nieprawidłowe wartości dla sekwencji protein")
        self.amino = protein
        self.length1 = len(self.amino)
        self.aminoFASTADict = self.FASTADict()

    def FASTADict(self):
        """
        Metoda zwraca słownik aminokwasów z odpowiednikami w formacie FASTA.

        Returns:
            Słownik z nazwami aminokwasów z odpowiednikami FASTA.
        """
        aminoDict = {
            "Ala": "A", "Asx": "B", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F", "Gly": "G",
            "His": "H", "Ile": "I", "Xle": "J", "Lys": "K", "Leu": "L", "Met": "M", "Asn": "N",
            "Pyl": "O", "Pro": "P", "Gln": "Q", "Arg": "R", "Ser": "S", "Thr": "T", "Sec": "U",
            "Val": "V", "Trp": "W", "Tyr": "Y", "Glx": "Z", "Ter": "*"
        }
        return aminoDict
    
    def Proteintext(self):
        """
        Metoda zwraca tekstową reprezentację sekwencji protein.

        Returns:
            Tekstowa reprezentacja sekwencji protein.
        """
        Text1 = f"> {self.identifier}\n {self.amino}"
        return Text1

    def ProteintextFASTA(self):
        """
        Metoda zwraca tekstową reprezentację sekwencji protein w formacie FASTA.

        Returns:
            Tekstowa reprezentacja sekwencji protein w formacie FASTA.
        """
        Text1 = f"> {self.identifier}\n"
        for elem in self.amino:
            Text1 += self.aminoFASTADict[elem]
        return Text1

    
    def mutate(self, position, value):
        """
        Metoda zamienia proteinę na zadanej pozycji sekwencji protein przez zastąpienie jej na określonej pozycji.

        Args:
            position: Pozycja, na której znajduje się wartość do zmiany.
            value: Nowa wartość, która zostanie wpisana w miejsce konkretnej pozycji.

        Raises:
            ValueError: W przypadku, gdy podany jest nieprawidłowa protein oraz, gdy liczba wykracza poza zakres listy protein.
        """
        if value not in self.valid_Chars:
            raise ValueError("Nieprawidłowa wartość dla sekwencji protein")
        if len(self.amino) < position:
            raise ValueError("Liczba wykraczająca poza zakres nici RNA")
        self.amino[position - 1] = value

    def findMotif(self, motif):
        """
        Metoda znajduje wszystkie pozycje, na których występuje zadany motyw w sekwencji protein.

        Args:
            motif: Motyw, który ma zostać znaleziony w sekwencji.

        Returns:
            Lista pozycji, gdzie motyw występuje w sekwencji.

        Raises:
            ValueError: W przypadku, gdy podana wartość jest nieprawidłowa dla sekwencji protein.
        """
        if motif not in self.valid_Chars:
            raise ValueError("Nieprawidłowe wartości dla sekwencji protein")
        listOfNucleotides = []
        for index, elem in enumerate(self.amino):
            if elem == motif:
                listOfNucleotides.append(index + 1)
        return listOfNucleotides
    

def main():
    """
    Główna funkcja programu.

    Prezentuje sposób działania programu na podstawie jednej nici.
    """
    DNAseqNew = DNASequence(1, "ATGAGGAAAGGGAAATAA")
    print("................................................")
    print("Wyświetlenie w formie tekstowej w formacie FASTA")
    print(DNAseqNew.DNAtext())
    DNAseqNew.mutate(6, 'A')
    print("................................................")
    print("Wyświetlenie sekwencji DNA po zmianie")
    print(DNAseqNew.DNAtext())
    newList = DNAseqNew.findMotif('A')
    print("................................................")
    print("Wyświetlenie pozycji zadanego motywu w sekwencji DNA w kierunku od 5' do 3'")
    print(newList)
    newDNAmat = DNAseqNew.Complement()
    print("................................................")
    print("Wyświetlenie komplementarnej matrycowej nici DNA w kierunku od 3' do 5' do nici "
          "kodującej DNA")
    print(newDNAmat)
    newRNA = DNAseqNew.transcribe(newDNAmat)
    print("................................................")
    print("Wyświetlenie komplementarnej nici RNA do nici matrycowej DNA ")
    print(newRNA.RNAtext())
    newRNA.mutate(6, 'U')
    print("................................................")
    print("Wyświetlenie sekwencji RNA po zmianie")
    print(newRNA.RNAtext())
    RNAlist = newRNA.findMotif("AAA")
    print("................................................")
    print("Wyświetlenie pozycji zadanego motywu w sekwencji RNA ")
    print(RNAlist)
    proteins = newRNA.transcribe()
    print("................................................")
    print("Wyświetlenie komplementarnej do RNA sekwencji protein ")
    print(proteins.Proteintext())
    proteins.mutate(2, "Cys")
    print("................................................")
    print("Wyświetlenie sekwencji protein po zmianie")
    print(proteins.Proteintext())
    proteinlist = proteins.findMotif("Lys")
    print("................................................")
    print("Wyświetlenie pozycji zadanego motywu w sekwencji protein ")
    print(proteinlist)
    print("................................................")
    print("Wyświetlenie sekwencji protein w formacie FASTA ")
    print(proteins.ProteintextFASTA())

if __name__ == "__main__":
    main()


import unittest

class DNASequenceTest(unittest.TestCase):

    def test_Constructor(self):
        validDNA = DNASequence(2, "AAAGGGCCTC")
        self.assertEqual("AAAGGGCCTC", validDNA.data, "Nieprawidłowe utworzenie obiektu klasy")
        self.assertEqual(2, validDNA.identifier, "Nieprawidłowe utworzenie obiektu klasy")

        with self.assertRaises(ValueError):
            DNASequence(1, "AAAGTTWRA")

    def test_DNAtext(self):
        DNAseq = DNASequence(1, "ATGAGGAAAGGGAAATAA")
        self.assertEqual("> 1\n ATGAGGAAAGGGAAATAA", DNAseq.DNAtext(), "Nieprawidłowe wyświetlenie nici DNA")

    def test_Mutate(self):
        DNAseq = DNASequence(1, "ATGAGGAAA")
        DNAseq.mutate(6, 'A')
        self.assertEqual("ATGAGAAAA", DNAseq.data, "Nieprawidłowe działanie metody mutate")

    def test_Mutate1(self):
        DNAseq = DNASequence(1, "ATGAGGAAA")
        with self.assertRaises(ValueError):
            DNAseq.mutate(6, 'W')

    def test_Mutate2(self):
        DNAseq = DNASequence(1, "ATGAGGAAA")
        with self.assertRaises(ValueError):
            DNAseq.mutate(24, 'G')

    def test_findMotif(self):
        DNAseq = DNASequence(1, "ATGAGGAAA")
        newList = DNAseq.findMotif('A')
        self.assertEqual("[1, 4, 7, 8, 9]", str(newList))

    def test_findMotif2(self):
        DNAseq = DNASequence(1, "ATGAGGAAA")
        with self.assertRaises(ValueError):
            DNAseq.findMotif('B')

    def test_Complement(self):
        DNAseq = DNASequence(1, "ATGAGGAAAGGGAAA")
        newDNAmat = DNAseq.Complement()
        self.assertEqual("TTTCCCTTTCCTCAT", newDNAmat, "Niepoprawnie utworzona nić komplementarna")

    def test_Transcribe(self):
        DNAseq = DNASequence(1, "ATGAGGAAG")
        newDNAmat = DNAseq.Complement()
        newRNA = DNAseq.transcribe(newDNAmat)
        self.assertEqual("AUGAGGAAG", newRNA.data, "Niepoprawna transkrypcja")


class RNASequenceTest(unittest.TestCase):

    def test_Constructor(self):
        validRNA = RNASequence(2, "AAAGGGCCU")
        self.assertEqual("AAAGGGCCU", validRNA.data, "Nieprawidłowe utworzenie obiektu klasy")
        self.assertEqual(2, validRNA.identifier, "Nieprawidłowe utworzenie obiektu klasy")

        with self.assertRaises(ValueError):
            RNASequence(1, "AAAGUUWR")

    def test_RNAtext(self):
        RNAseq = RNASequence(1, "AGAGGAAAGGGAAAA")
        self.assertEqual("> 1\n AGAGGAAAGGGAAAA", RNAseq.RNAtext(), "Nieprawidłowe wyświetlenie nici RNA")

    def test_Mutate(self):
        RNAseq = RNASequence(1, "AUGAGGAAA")
        RNAseq.mutate(6, 'A')
        self.assertEqual("AUGAGAAAA", RNAseq.data, "Nieprawidłowe działanie metody mutate")

    def test_Mutate1(self):
        RNAseq = RNASequence(1, "AUGAGGAAA")
        with self.assertRaises(ValueError):
            RNAseq.mutate(6, 'W')

    def test_Mutate2(self):
        RNAseq = RNASequence(1, "AUGAGGAAA")
        with self.assertRaises(ValueError):
            RNAseq.mutate(24, 'G')

    def test_findMotif(self):
        RNAseq = RNASequence(1, "AUGAGGAAA")
        newList = RNAseq.findMotif("AUG")
        self.assertEqual("[1]", str(newList), "Nieprawidłowo obliczona pozycja")

    def test_findMotif2(self):
        RNAseq = RNASequence(1, "AUGAGGAAA")
        with self.assertRaises(ValueError):
            RNAseq.findMotif("B")

    def test_findMotif3(self):
        RNAseq = RNASequence(1, "AUGAGGAAAG")
        with self.assertRaises(ValueError):
            RNAseq.findMotif("AAA")

    def test_Transcribe(self):
        RNAseq = RNASequence(1, "AUGAGGAAGUAA")
        newProtein = RNAseq.transcribe()
        self.assertEqual("['Met', 'Arg', 'Lys', 'Ter']", str(newProtein.amino), "Niepoprawna transkrypcja")

    def test_Transcribe1(self):
        RNAseq = RNASequence(1, "AGGAAGUAA")
        with self.assertRaises(ValueError):
            RNAseq.transcribe()

    def test_Transcribe2(self):
        RNAseq = RNASequence(1, "AUGAGGAAG")
        with self.assertRaises(ValueError):
            RNAseq.transcribe()


class ProteinSequenceTest(unittest.TestCase):

    def test_Constructor(self):
        validprotein = ProteinSequence(2, ["Met", "Arg", "Lys", "Ter"])
        self.assertEqual(["Met", "Arg", "Lys", "Ter"], validprotein.amino, "Nieprawidłowe utworzenie obiektu klasy")
        self.assertEqual(2, validprotein.identifier, "Nieprawidłowe utworzenie obiektu klasy")

        with self.assertRaises(ValueError):
            ProteinSequence(1, ["Met", "Arg", "Asa", "Lys", "Ter"])

    def test_proteintext(self):
        proteinseq = ProteinSequence(1, ["Met", "Arg", "Lys", "Ter"])
        self.assertEqual("> 1\n ['Met', 'Arg', 'Lys', 'Ter']", proteinseq.Proteintext(), "Nieprawidłowe wyświetlenie protein")

    def test_proteintextFASTA(self):
        proteinseq = ProteinSequence(1, ["Met", "Arg", "Lys", "Ter"])
        self.assertEqual("> 1\nMRK*", proteinseq.ProteintextFASTA(), "Nieprawidłowe wyświetlenie protein w formacie FASTA")

    def test_Mutate(self):
        proteinseq = ProteinSequence(1, ["Met", "Arg", "Lys", "Ter"])
        proteinseq.mutate(2, "Lys")
        self.assertEqual(["Met", "Lys", "Lys", "Ter"], proteinseq.amino, "Nieprawidłowe działanie metody mutate")

    def test_Mutate1(self):
        proteinseq = ProteinSequence(1, ["Met", "Arg", "Lys", "Ter"])
        with self.assertRaises(ValueError):
            proteinseq.mutate(3, "Was")

    def test_Mutate2(self):
        proteinseq = ProteinSequence(1, ["Met", "Arg", "Lys", "Ter"])
        with self.assertRaises(ValueError):
            proteinseq.mutate(24, "Lys")

    def test_findMotif(self):
        proteinseq = ProteinSequence(1, ["Met", "Arg", "Lys", "Ter"])
        newList = proteinseq.findMotif("Arg")
        self.assertEqual("[2]", str(newList), "Nieprawidłowo obliczona pozycja")

    def test_findMotif2(self):
        proteinseq = ProteinSequence(1, ["Met", "Arg", "Lys", "Ter"])
        with self.assertRaises(ValueError):
            proteinseq.findMotif("Byh")


if __name__ == "__main__":
    unittest.main()

