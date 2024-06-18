import math

class Wielomian:
    """
    Klasa reprezentująca wielomian dowolnego stopnia.

        Attributes:
            listOfFactors: Lista zawierająca współczynniki wielomianu.

        Args:
            listOfValue: Lista zawierająca współczynniki wielomianu.

        Raises:
            ValueError: W przypadku, gdy lista współczynników jest pusta.
    """
    def __init__(self, listOfValue):
        if not listOfValue:
            raise ValueError("Próba niepoprawnego utworzenia wielomianu")
        
        while listOfValue[0] == 0.0:
            listOfValue.pop(0)
        
        self.listOfFactors = listOfValue

    def StopienWielomianu(self):
        """
        Metoda zwraca stopień wielomianu.

        Returns:
            Stopień wielomianu.
        """
        return len(self.listOfFactors) - 1

    def TekstWielomian(self):
        """
        Metoda generuje tekstową reprezentację wielomianu.

        Returns:
            Tekstowa reprezentacja wielomianu.
        """
        temp_Factor = self.StopienWielomianu()
        resultOfPolynomial = "W(x)= "
        firstElem = True
        
        for index, elem in enumerate(self.listOfFactors):
            if index == len(self.listOfFactors) - 1:
                if elem == 0.0:
                    break
                if elem > 0.0:
                    resultOfPolynomial += '+' + str(elem)
                else:
                    resultOfPolynomial += str(elem)
                break
            
            if elem != 0.0:
                if elem > 0.0:
                    if firstElem:
                        firstElem = False
                    else:
                        resultOfPolynomial += '+'
                if firstElem:
                    firstElem = False
                
                resultOfPolynomial += str(elem)
                if index == len(self.listOfFactors) - 2:
                    resultOfPolynomial += 'x'
                else:
                    resultOfPolynomial += f"x^{temp_Factor}"
            temp_Factor -= 1

        return resultOfPolynomial

    def __call__(self, argument1):
        """
        Przeciązenie operatora wywołania oblicza wartość wielomianu dla danego argumentu.

        Args:
            argument1: Argument dla którego obliczana jest wartość wielomianu.

        Returns:
            Wynik obliczenia wielomianu dla danego argumentu.
        """
        result = 0.0
        tempFactor1 = self.StopienWielomianu()
        
        for elem in self.listOfFactors:
            result += argument1 ** tempFactor1 * elem
            tempFactor1 -= 1
        
        return result

    def __add__(self, polynomial_2):
        """
        Przeciązony operator dodawania wielomianów.

        Args:
            polynomial_2: Drugi wielomian do dodania.

        Returns:
            Nowy wielomian będący wynikiem dodawania.
        """
        listW = []
        degreeW1 = self.StopienWielomianu()
        degreeW2 = polynomial_2.StopienWielomianu()
        difference = abs(degreeW1 - degreeW2)
        
        if degreeW1 != degreeW2:
            if degreeW1 > degreeW2:
                listW = self.listOfFactors[:difference]
                for i in range(degreeW2 + 1):
                    listW.append(self.listOfFactors[difference + i] + polynomial_2.listOfFactors[i])
            else:
                listW = polynomial_2.listOfFactors[:difference]
                for i in range(degreeW1 + 1):
                    listW.append(polynomial_2.listOfFactors[difference + i] + self.listOfFactors[i])
        else:
            for i in range(degreeW1 + 1):
                listW.append(polynomial_2.listOfFactors[i] + self.listOfFactors[i])

        polynomial3 = Wielomian(listW)
        return polynomial3

    def __iadd__(self, polynomial_2):
        """
        Przeciązony operator złozony dodawania wielomianów.

        Args:
            polynomial_2: Drugi wielomian do dodania.

        Returns:
            Referencja do obiektu bazowego po wykonaniu dodawania.
        """
        listW = []
        degreeW1 = self.StopienWielomianu()
        degreeW2 = polynomial_2.StopienWielomianu()
        difference = abs(degreeW1 - degreeW2)
        
        if degreeW1 != degreeW2:
            if degreeW1 > degreeW2:
                listW = self.listOfFactors[:difference]
                for i in range(degreeW2 + 1):
                    listW.append(self.listOfFactors[difference + i] + polynomial_2.listOfFactors[i])
            else:
                listW = polynomial_2.listOfFactors[:difference]
                for i in range(degreeW1 + 1):
                    listW.append(polynomial_2.listOfFactors[difference + i] + self.listOfFactors[i])
        else:
            for i in range(degreeW1 + 1):
                listW.append(polynomial_2.listOfFactors[i] + self.listOfFactors[i])

        self.listOfFactors = listW
        return self

    def __sub__(self, polynomial_2):
        """
        Przeciązony operator odejmowania wielomianów.

        Args:
            polynomial_2: Drugi wielomian do odjęcia.

        Returns:
            Nowy wielomian będący wynikiem odejmowania.
        """
        listW = []
        degreeW1 = self.StopienWielomianu()
        degreeW2 = polynomial_2.StopienWielomianu()
        difference = abs(degreeW1 - degreeW2)
        
        if degreeW1 != degreeW2:
            if degreeW1 > degreeW2:
                listW = self.listOfFactors[:difference]
                for i in range(degreeW2 + 1):
                    listW.append(self.listOfFactors[difference + i] - polynomial_2.listOfFactors[i])
            else:
                listW = polynomial_2.listOfFactors[:difference]
                for i in range(len(listW)):
                    listW[i] *= -1.0
                for i in range(degreeW1 + 1):
                    listW.append(self.listOfFactors[i] - polynomial_2.listOfFactors[difference + i])
        else:
            for i in range(degreeW1 + 1):
                listW.append(self.listOfFactors[i] - polynomial_2.listOfFactors[i])

        polynomial3 = Wielomian(listW)
        return polynomial3

    def __isub__(self, polynomial_2):
        """
        Przeciązony operator złozony odejmowania wielomianów.

        Args:
            polynomial_2: Drugi wielomian do odjęcia.

        Returns:
            Referencja do obiektu bazowego po odejmowaniu.
        """
        listW = []
        degreeW1 = self.StopienWielomianu()
        degreeW2 = polynomial_2.StopienWielomianu()
        difference = abs(degreeW1 - degreeW2)
        
        if degreeW1 != degreeW2:
            if degreeW1 > degreeW2:
                listW = self.listOfFactors[:difference]
                for i in range(degreeW2 + 1):
                    listW.append(self.listOfFactors[difference + i] - polynomial_2.listOfFactors[i])
            else:
                listW = polynomial_2.listOfFactors[:difference]
                for i in range(len(listW)):
                    listW[i] *= -1.0
                for i in range(degreeW1 + 1):
                    listW.append(self.listOfFactors[i] - polynomial_2.listOfFactors[difference + i])
        else:
            for i in range(degreeW1 + 1):
                listW.append(self.listOfFactors[i] - polynomial_2.listOfFactors[i])

        self.listOfFactors = listW
        return self

    def __mul__(self, polynomial_2):
        """
        Przeciązony operator mnożenia wielomianów.

        Args:
            polynomial_2: Drugi wielomian do przemnożenia.

        Returns:
            Nowy wielomian będący wynikiem mnożenia.
        """
        x = len(self.listOfFactors) + len(polynomial_2.listOfFactors) - 1
        listW = [0.0] * x
        
        for index1, elem in enumerate(self.listOfFactors):
            for index2, elem2 in enumerate(polynomial_2.listOfFactors):
                listW[index1 + index2] += self.listOfFactors[index1] * polynomial_2.listOfFactors[index2]

        return Wielomian(listW)

    def __imul__(self, polynomial_2):
        """
        Przeciązony operator złozony mnożenia wielomianów.

        Args:
            polynomial_2: Drugi wielomian do pomnożenia.

        Returns:
            Referencja do obiektu bazowego po wykonaniu mnożenia.
        """
        x = len(self.listOfFactors) + len(polynomial_2.listOfFactors) - 1
        listW = [0.0] * x
        
        for index1, elem in enumerate(self.listOfFactors):
            for index2, elem2 in enumerate(polynomial_2.listOfFactors):
                listW[index1 + index2] += self.listOfFactors[index1] * polynomial_2.listOfFactors[index2]

        self.listOfFactors = listW
        return self


def main():
    """
    Funkcja main służąca do demonstracji operacji na wielomianach.

    Funkcja tworzy wielomian na podstawie listy wartości.
    Funkcja oblicza stopień wielomianu oraz wyświetla jego tekstową reprezentację.
    Funkcja oblicza wartość wielomianu dla określonej wartości.
    Funkcja tworzy drugi wielomian.
    Funkcja dodaje dwa wielomiany i wyświetla wynik.
    Funkcja odejmuje drugi wielomian od pierwszego i wyświetla wynik.
    Funkcja mnoży dwa wielomiany i wyświetla wynik.
    """
    listOfValue = [-4.0, 2.0]
    newPolynomial = Wielomian(listOfValue)
    degreeOfPolynomial = newPolynomial.StopienWielomianu()
    print("................................................")
    print("Stopień wielomianu:", degreeOfPolynomial)
    print("................................................")
    print("Wielomian pierwszy: ")
    textPolynomial1 = newPolynomial.TekstWielomian()
    print(textPolynomial1)
    resultOfPolynomial = newPolynomial(9.0)
    print("................................................")
    print("Wynik wielomianu:", resultOfPolynomial)

    Polynomial_2 = Wielomian([0.0, 0.5, 2.0])
    print("................................................")
    print("Wielomian drugi: ")
    print(Polynomial_2.TekstWielomian())

    addResult = newPolynomial + Polynomial_2
    print("................................................")
    print("Wynik doodawania:", newPolynomial.TekstWielomian(), "+", Polynomial_2.TekstWielomian())
    print(addResult.TekstWielomian())

    minusResult = newPolynomial - Polynomial_2
    print("................................................")
    print("Wynik odejmowania:", newPolynomial.TekstWielomian(), "-", Polynomial_2.TekstWielomian())
    print(minusResult.TekstWielomian())

    multiplicationResult = newPolynomial * Polynomial_2
    print("................................................")
    print("Wynik mnożenia:", newPolynomial.TekstWielomian(), "*", Polynomial_2.TekstWielomian())
    print(multiplicationResult.TekstWielomian())


if __name__ == "__main__":
    main()


import unittest

class WielomianTest(unittest.TestCase):

    def testConstructor(self):
        validPolynomial = Wielomian([1.0, 2.0, 3.0])
        self.assertEqual([1.0, 2.0, 3.0], validPolynomial.listOfFactors, "Niepoprawnie utworzone współczynniki wielomianu")

        zeroPolynomial = Wielomian([0.0, 0.0, 2.0])
        self.assertEqual([2.0], zeroPolynomial.listOfFactors, "Niepoprawnie utworzone współczynniki wielomianu")

        with self.assertRaises(ValueError):
            Wielomian([])

    def test_polynomialDegree(self):
        wielomian = Wielomian([2.0, 3.0, 4.0])
        self.assertEqual(2, wielomian.StopienWielomianu(), "Niepoprawnie obliczony stopień wielomianu")

    def test_polynomialDegree1(self):
        wielomian = Wielomian([2.0, 3.0, 0.0, 4.0])
        self.assertEqual(3, wielomian.StopienWielomianu(), "Niepoprawnie obliczony stopień wielomianu")

    def test_polynomialDegree2(self):
        wielomian = Wielomian([0.0, 0.0, 2.0, 3.0, 0.0, 4.0])
        self.assertEqual(3, wielomian.StopienWielomianu(), "Niepoprawnie obliczony stopień wielomianu")

    def test_polynomialText(self):
        listOfValue = [-4.0, 2.0, 0.0, 0.0]
        wielomian = Wielomian(listOfValue)
        tekstWielomianu = wielomian.TekstWielomian()
        self.assertEqual("W(x)= -4.0x^3+2.0x^2", tekstWielomianu, "Zły zapis wielomianu")

    def test_resultOfPolynomailInvoke(self):
        listOfValue = [-4.0, 2.0, 0.0, 0.0]
        wielomian = Wielomian(listOfValue)
        result = wielomian(2.0)
        self.assertEqual(-24.0, result, "Zły wynik wielomianu")

    def test_polynomialAdd(self):
        wielomian1 = Wielomian([-4.0, 2.0, 0.0, 0.0])
        wielomian2 = Wielomian([0.0, 0.5, 2.0])
        wynikDodawania = wielomian1 + wielomian2
        self.assertEqual("W(x)= -4.0x^3+2.0x^2+0.5x+2.0", wynikDodawania.TekstWielomian(), "Niepoprawny wynik dodawania")

    def test_polynomialMinus(self):
        wielomian1 = Wielomian([-4.0, 2.0, 0.0, 0.0])
        wielomian2 = Wielomian([0.0, 0.5, 2.0])
        wynikOdejmowania = wielomian1 - wielomian2
        self.assertEqual("W(x)= -4.0x^3+2.0x^2-0.5x-2.0", wynikOdejmowania.TekstWielomian(), "Niepoprawny wynik odejmowania")

    def test_polynomialMultiplication(self):
        wielomian1 = Wielomian([-4.0, 2.0, 0.0, 0.0])
        wielomian2 = Wielomian([0.0, 0.5, 2.0])
        wynikMnozenia = wielomian1 * wielomian2
        self.assertEqual("W(x)= -2.0x^4-7.0x^3+4.0x^2", wynikMnozenia.TekstWielomian(), "Niepoprawny wynik mnożenia")


if __name__ == "__main__":
    unittest.main()
