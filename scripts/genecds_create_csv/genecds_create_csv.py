import subprocess


def launch_gene_cds_field_distribution():
    print("Launch genecds_field_distribution.py\n")
    subprocess.Popen(r"python genecds_field_distribution.py "
                     r"-input Mamastrovirus_1_complete_genome_records.gb "
                     r"-oname Mamastrovirus_1_gene_product_output.txt").communicate()
    print("\nOutput saved to:\nMamastrovirus_1_gene_product_output.txt\n")


def gene_cds_field_distribution_output_handle():

    print("Creating csv file...\n")
    with open(r"Mamastrovirus_1_gene_product_output.txt", "r") as output:
        with open(r"Mamastrovirus_1_CDS_gene_product.csv", "w") as writer:
            for line in output.readlines()[4:]:
                if (l := line.strip().lower()) == "orf1" or l == "orf1a" \
                        or "unknown" in l or "protease" in l or "NS" in line.strip():
                    print(f"{line.strip()},1A", file=writer)
                elif "polymerase" in l or l == "orf1b":
                    print(f"{line.strip()},1B", file=writer)
                elif "Non" in line.strip():
                    print(f"{line.strip()},1AB", file=writer)
                elif "non" in line:
                    print(f"{line.strip()},1A\n"
                          f"{line.strip()},1B\n"
                          f"{line.strip()},1AB", file=writer)
                else:
                    print(f"{line.strip()},2", file=writer)
    print("Done")


launch_gene_cds_field_distribution()
gene_cds_field_distribution_output_handle()
