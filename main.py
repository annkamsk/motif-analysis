import csv
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio import motifs

from scipy import stats


PROTEIN_FRAGMENTS_PATH = "data/protein_fragments.fa"
ECOLI_PATH = "data/genes_e_coli_new.fa"
ECOLI_TRANS = "data/protein_e_coli.fa"
ECOLI_OUT = "data/ecoli.xml"
BLASTP = "~/Programs/ncbi-blast-2.10.0+/bin/blastp"
FOUND_GENES = "data/genes.fa"
PROMS = "data/proms_e_coli_fixed.fa"
PROMS_A = "data/promsA.fa"
PROMS_B = "data/promsB.fa"
MEMEA = "data/memeA.xml"
MEMEB = "data/memeB.xml"
TASK1 = "data/task1.csv"
TASK2_A = "data/task2_A.pfm"
TASK2_B = "data/task2_B.pfm"
TASK3 = "data/task3.csv"


def translate_ecoli():
    ecoli = []
    for record in SeqIO.parse(ECOLI_PATH, "fasta"):
        ecoli.append(
            SeqRecord(
                record.seq.translate(to_stop=True),
                id=record.id,
                description=record.description,
            )
        )

    with open(ECOLI_TRANS, "w+") as output:
        SeqIO.write(ecoli, output, "fasta")


def hit_db():
    blastx_cline = NcbiblastxCommandline(
        cmd=BLASTP,
        query=PROTEIN_FRAGMENTS_PATH,
        db=ECOLI_TRANS,
        evalue=0.001,
        outfmt=5,
        out=ECOLI_OUT,
    )
    blastx_cline()


def parse_blast():
    with open(ECOLI_OUT) as handle:
        blast_records = list(NCBIXML.parse(handle))

    matches = []
    for record in blast_records:
        best = max(
            zip(record.alignments, record.descriptions), key=lambda p: p[1].score
        )[0]
        for hsp in best.hsps:
            matches.append(
                {
                    "input_id": record.query,
                    "gene_id": best.title.split()[0],
                    "e-value": hsp.expect,
                }
            )
    return matches


def save_result(data, path):
    with open(path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        for d in data:
            if isinstance(d, type({})):
                d = [*d.values()]
            writer.writerow(d)


def read_data(path):
    proteins = []
    for record in SeqIO.parse(path, "fasta"):
        proteins.append(
            SeqRecord(record.seq, id=record.id, description=record.description)
        )
    return proteins


def get_hits(m, group):
    hits = 0
    total = 0
    # obtain a number of positions in these promoter sequences that have a log-odds score higher than 0
    for prom in group:
        total += len(prom.seq) - m.length
        pssm = m.counts.normalize().log_odds()
        hits += len(list(pssm.search(prom.seq, threshold=0.0)))
    return hits, total


def main():
    """
    TASK 1:
    """
    translate_ecoli()

    # create local BLAST db
    os.system(
        'makeblastdb -in data/protein_e_coli.fa -parse_seqids -blastdb_version 5 -taxid 1 -title "ecoli" -dbtype prot'
    )

    hit_db()
    matches = parse_blast()
    save_result(matches, TASK1)

    """
        TASK 2:
    """
    # gene to group
    proms = read_data(PROMS)
    for match in matches:
        match["prom"] = [prom for prom in proms if prom.id == match["gene_id"]][0]

    # separate group A and B and pack into dict for unique values
    promsA = [
        *dict(
            (m["prom"].id, m["prom"])
            for m in matches
            if m["input_id"].startswith("groupA")
        ).values()
    ]
    promsB = [
        *dict(
            (m["prom"].id, m["prom"])
            for m in matches
            if m["input_id"].startswith("groupB")
        ).values()
    ]

    # write groups of proms to file, so they can be sent to MEME
    with open(PROMS_A, "w+") as output:
        SeqIO.write(promsA, output, "fasta")
    with open(PROMS_B, "w+") as output:
        SeqIO.write(promsB, output, "fasta")

    # parse results from MEME
    with open(MEMEA, "r") as handle:
        recordA = motifs.parse(handle, "meme")
    with open(MEMEB, "r") as handle:
        recordB = motifs.parse(handle, "meme")

    # save result
    with open(TASK2_A, "w+") as output:
        for m in recordA:
            output.write(m.format("pfm"))
            output.write("\n")

    with open(TASK2_B, "w+") as output:
        for m in recordB:
            output.write(m.format("pfm"))
            output.write("\n")

    """
        TASK 3:    
    """
    result = []
    for m in recordA:
        hitsA, totalA = get_hits(m, promsA)
        hitsB, totalB = get_hits(m, promsB)
        p = stats.binom_test(hitsA, n=totalA, p=hitsB / totalB)
        result.append(["A", hitsA, p])

    for m in recordB:
        hitsA, totalA = get_hits(m, promsA)
        hitsB, totalB = get_hits(m, promsB)
        p = stats.binom_test(hitsB, n=totalB, p=hitsA / totalA)
        result.append(["B", hitsB, p])

    save_result(result, TASK3)


if __name__ == "__main__":
    main()
