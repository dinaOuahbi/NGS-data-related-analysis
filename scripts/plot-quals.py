import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)
plt.xlabel("La confiance dans l'exactitude du variant appélé (QUAL)")
plt.ylabel('Frequence')

plt.savefig(snakemake.output[0])
