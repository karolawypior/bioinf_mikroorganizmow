from Bio import SeqIO
import statistics

def zlicz_geny(plik_gbff):
    licznik_znane = 0
    licznik_nieznane = 0
    lista_dlugosci = []

    genbank = SeqIO.parse(plik_gbff, "genbank")

    for rekord in genbank:
        cechy = rekord.features
        for cecha in cechy:
            if cecha.type != "CDS":
                continue

            start = int(cecha.location.start)
            koniec = int(cecha.location.end)
            dlugosc = koniec - start
            lista_dlugosci.append(dlugosc)

            produkty = cecha.qualifiers.get("product", [""])
            produkt = produkty[0].lower() if produkty else ""

            if produkt == "" or "hypothetical" in produkt or "unknown" in produkt:
                licznik_nieznane += 1
            else:
                licznik_znane += 1

    mediana = statistics.median(lista_dlugosci) if lista_dlugosci else 0
    return licznik_znane, licznik_nieznane, mediana


def pokaz_wyniki(sciezka, opis):
    znane, nieznane, mediana = zlicz_geny(sciezka)
    print(f"{opis}:")
    print(f"  Geny o znanej funkcji: {znane}")
    print(f"  Geny o nieznanej funkcji: {nieznane}")
    print(f"  Mediana długości genu: {mediana}\n")


def main():
    plik_light = "results_bakta_light/assembly.gbff"
    plik_full = "bakta_full/assembly.gbff"

    pokaz_wyniki(plik_light, "Bakta light")
    pokaz_wyniki(plik_full, "Bakta full")
    
if __name__ == "__main__":
    main()

