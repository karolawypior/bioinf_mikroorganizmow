#!/usr/bin/env python3

import gzip
import sys
import statistics
import matplotlib.pyplot as plt


def l_odczytow_fastq_gz(fastq_gz):
    """Liczy liczbę odczytów w FASTQ.gz"""
    with gzip.open(fastq_gz, "rt") as f:
        il_linii = sum(1 for _ in f)
    return il_linii // 4


def statystyki_nanopore(fastq_gz):
    """Zbiera statystyki Nanopore"""
    dlugosc = []
    l_nukleotydow = 0
    gc_licznik = 0

    with gzip.open(fastq_gz, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # linia sekwencji
                seq = line.strip()
                l = len(seq)
                dlugosc.append(l)
                l_nukleotydow += l
                gc_licznik += seq.count("G") + seq.count("C")

    return dlugosc, l_nukleotydow, gc_licznik


def plot_histogram(dlugosc, output="nanopore_read_length_histogram.png"):
    plt.hist(dlugosc, bins=50)
    plt.xlabel("Długość odczytu (bp)")
    plt.ylabel("Liczba odczytów")
    plt.title("Rozkład długości odczytów Nanopore")
    plt.savefig(output, dpi=300)
    plt.close()


def main():
    if len(sys.argv) != 4:
        print("Użycie: python zad1_lab1.py IL_1.fastq.gz IL_2.fastq.gz NP.fastq.gz")
        sys.exit(1)

    IL_1 = sys.argv[1]
    IL_2 = sys.argv[2]
    nanopore = sys.argv[3]

    # --- Illumina ---
    odczyty_IL1 = l_odczytow_fastq_gz(IL_1)
    odczyty_IL2 = l_odczytow_fastq_gz(IL_2)
    total_illumina = odczyty_IL1 + odczyty_IL2

    print("ILLUMINA")
    print(f"  IL1 odczyty: {odczyty_IL1}")
    print(f"  IL2 odczyty: {odczyty_IL2}")
    print(f"  Całkowita liczba odczytów Illumina: {total_illumina}")
    print()

    # --- Nanopore ---
    dlugosc,l_nukleotydow, gc_licznik = statystyki_nanopore(nanopore)

    print("NANOPORE")
    print(f"  Całkowita liczba nukleotydów: {l_nukleotydow}")
    print(f"  Średnia długość odczytu: {statistics.mean(dlugosc):.2f}")
    print(f"  Mediana długości odczytu: {statistics.median(dlugosc):.2f}")

    if l_nukleotydow > 0:
        gc_procent = (gc_licznik / l_nukleotydow) * 100
        print(f"  %GC: {gc_procent:.2f}%")

    plot_histogram(dlugosc)
    print("\nHistogram zapisano jako: nanopore_read_length_histogram.png")


if __name__ == "__main__":
    main()
