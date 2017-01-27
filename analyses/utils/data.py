import pandas as pd
from os import path, getcwd, environ
import numpy as np

from cohorts import Sample, Patient, Cohort, DataFrameLoader
from cohorts.variant_stats import variant_stats_from_variant

# This repo's "data" directory
REPO_DATA_DIR = path.join(path.dirname(path.dirname(getcwd())), "data")

def get_dir(env_var, default_dir):
    dir = environ.get(env_var, None)
    if dir is None:
        return default_dir
    return dir

# The directory where we want to store cached mutations, neoantigens, etc,
DATA_DIR = "/demeter/scratch/datasets/mskcc/bladder/"
CACHE_DATA_DIR = get_dir("CACHE_DATA_DIR", path.join(REPO_DATA_DIR, "cache"))
PAGEANT_COVERAGE_DATA_DIR = get_dir("PAGEANT_COVERAGE_DATA_DIR", path.join(REPO_DATA_DIR, "pageant_ensembl_coverage_output"))
VCF_DATA_DIR = get_dir("VCF_DATA_DIR", path.join(DATA_DIR, "vcfs-indelrealigned-bqsr"))
BAM_DNA_DATA_DIR = get_dir("BAM_DNA_DATA_DIR", path.join(DATA_DIR, "bams"))
BAM_RNA_DATA_DIR = get_dir("BAM_RNA_DATA_DIR", path.join(DATA_DIR, "rna-aligned/star-aligned-bams"))
KALLISTO_DATA_DIR = get_dir("KALLISTO_DATA_DIR", path.join(DATA_DIR, "rna-aligned/kallisto"))
CUFFLINKS_DATA_DIR = get_dir("CUFFLINKS_DATA_DIR", path.join(DATA_DIR, "rna-aligned/star-cufflinks"))
TCR_DATA_DIR = get_dir("TCR_DATA_DIR", path.join(DATA_DIR, "tcr"))
POLYPHEN_DATA_DUMP = get_dir("POLYPHEN_DATA_DUMP", path.join(DATA_DIR, "polyphen-2.2.2-whess-2011_12.sqlite"))

def mutect_snv_file_format_func(patient_id, normal_bam_id, tumor_bam_id):
    return "Mutect-MSKCC-Bladder-%d-normal=%s.bam-indelrealigned-default-default-bqsr-default-default.bam-tumor=%s.bam-indelrealigned-default-default-bqsr-default-default.bam-merged.vcf" % (
        int(patient_id), normal_bam_id, tumor_bam_id)

def strelka_snv_file_format_func(patient_id, normal_bam_id, tumor_bam_id):
    return "Strelka-MSKCC-Bladder-%d-normal=%s.bam-indelrealigned-default-default-bqsr-default-default.bam-tumor=%s.bam-indelrealigned-default-default-bqsr-default-default.bamstrelka_output/results/passed.somatic.snvs.vcf" % (
        int(patient_id), normal_bam_id, tumor_bam_id)

def strelka_indel_file_format_func(patient_id, normal_bam_id, tumor_bam_id):
    return "Strelka-MSKCC-Bladder-%d-normal=%s.bam-indelrealigned-default-default-bqsr-default-default.bam-tumor=%s.bam-indelrealigned-default-default-bqsr-default-default.bamstrelka_output/results/passed.somatic.indels.vcf" % (
        int(patient_id), normal_bam_id, tumor_bam_id)

def mutect_without_bqsr_snv_file_format_func(patient_id, normal_bam_id, tumor_bam_id):
    return "Mutect-MSKCC-Bladder-%d-normal=%s.bam-tumor=%s.bam-merged.vcf" % (
        int(patient_id), normal_bam_id, tumor_bam_id)

def strelka_without_bqsr_snv_file_format_func(patient_id, normal_bam_id, tumor_bam_id):
    return "Strelka-MSKCC-Bladder-%d-normal=%s.bam-tumor=%s.bamstrelka_output/results/passed.somatic.snvs.vcf"  % (
        int(patient_id), normal_bam_id, tumor_bam_id)

def rizvi_filter(filterable_variant):
    somatic_stats = variant_stats_from_variant(filterable_variant.variant,
                                               filterable_variant.variant_metadata)
    return (
        somatic_stats.tumor_stats.depth >= 7 and
        somatic_stats.normal_stats.depth >= 7 and
        somatic_stats.tumor_stats.variant_allele_frequency > 0.1 and
        (1 - somatic_stats.normal_stats.variant_allele_frequency) > 0.97
    )

def init_cohort(**kwargs):
    return BladderData(**kwargs).load_cohort()

class BladderData(object):
    """
    Represents the bladder cohort.

    Parameters
    ----------
    exclude_patient_ids : set
        Exclude a particular collection of patient IDs from the cohort.
    only_patients_with_bams : bool
        If True, only include patients that have sequencing data.
    join_with : list
        Which dataframes to join on by default.
    join_how : list
        How to join those dataframes by default.
    """
    def __init__(self,
                 exclude_patient_ids={"4072"}, # 4072 had some failed sequencing
                 only_patients_with_bams=True,
                 clinical_csv_filename="clinical_updated.csv",
                 join_with=[],
                 join_how="inner",
                 filter_fn=rizvi_filter,
                 normalized_per_mb=True,
                 min_coverage_normal_depth=7,
                 min_coverage_tumor_depth=7,
                 responder_pfs_equals_os=False,
                 repo_data_dir=REPO_DATA_DIR,
                 cache_data_dir=CACHE_DATA_DIR,
                 vcf_data_dir=VCF_DATA_DIR,
                 bam_dna_data_dir=BAM_DNA_DATA_DIR,
                 bam_rna_data_dir=BAM_RNA_DATA_DIR,
                 kallisto_data_dir=KALLISTO_DATA_DIR,
                 cufflinks_data_dir=CUFFLINKS_DATA_DIR,
                 tcr_data_dir=TCR_DATA_DIR,
                 polyphen_data_dump=POLYPHEN_DATA_DUMP,
                 pageant_coverage_data_dir=PAGEANT_COVERAGE_DATA_DIR,
                 without_bqsr=False,
                 benefit_plot_name="DCB",
                 benefit_os_plot_name="DCB-OS",
                 hazard_plot_name="Progression or Mortality",
                 hazard_os_plot_name="Mortality",
                 tcr_peripheral_a_clonality_plot_name="Pre-treatment Peripheral TCR Clonality",
                 tcr_peripheral_a_clonality_short_plot_name="Pre-treatment Clonality",
                 tcr_tumor_clonality_plot_name="TIL Clonality",
                 til_fraction_plot_name="TIL Proportion",
                 **kwargs):
        self.exclude_patient_ids = exclude_patient_ids
        self.only_patients_with_bams = only_patients_with_bams
        self.clinical_csv_filename = clinical_csv_filename
        self.join_with = join_with
        self.join_how = join_how
        self.filter_fn = filter_fn
        self.normalized_per_mb = normalized_per_mb
        self.min_coverage_normal_depth = min_coverage_normal_depth
        self.min_coverage_tumor_depth = min_coverage_tumor_depth
        self.responder_pfs_equals_os = responder_pfs_equals_os
        self.repo_data_dir = repo_data_dir
        self.cache_data_dir = cache_data_dir
        self.vcf_data_dir = vcf_data_dir
        self.bam_dna_data_dir = bam_dna_data_dir
        self.bam_rna_data_dir = bam_rna_data_dir
        self.kallisto_data_dir = kallisto_data_dir
        self.cufflinks_data_dir = cufflinks_data_dir
        self.tcr_data_dir = tcr_data_dir
        self.polyphen_data_dump = polyphen_data_dump
        self.pageant_coverage_data_dir = pageant_coverage_data_dir
        self.without_bqsr = without_bqsr
        self.benefit_plot_name = benefit_plot_name
        self.benefit_os_plot_name = benefit_os_plot_name
        self.hazard_plot_name = hazard_plot_name
        self.hazard_os_plot_name = hazard_os_plot_name
        self.tcr_peripheral_a_clonality_plot_name = tcr_peripheral_a_clonality_plot_name
        self.tcr_peripheral_a_clonality_short_plot_name = tcr_peripheral_a_clonality_short_plot_name
        self.tcr_tumor_clonality_plot_name = tcr_tumor_clonality_plot_name
        self.til_fraction_plot_name = til_fraction_plot_name

        self.init_kwargs = kwargs

        if without_bqsr:
            vcf_data_dir_override = path.join(get_dir("DATA_DIR", DATA_DIR), "vcfs")
            print("Overriding VCF dir (%s) to use without-BQSR VCFs: %s" % (self.vcf_data_dir, vcf_data_dir_override))
            self.vcf_data_dir = vcf_data_dir_override

        def plot_benefit_os(self, **kwargs):
            no_benefit_os_plot_name = "No %s" % benefit_os_plot_name
            return self.plot_benefit(
                benefit_col="is_benefit_os",
                boolean_value_map={True: benefit_os_plot_name, False: no_benefit_os_plot_name},
                order=[no_benefit_os_plot_name, benefit_os_plot_name],
                **kwargs)
        Cohort.plot_benefit_os = plot_benefit_os

    def load_cohort(self, join_with=None, join_how="outer"):
        df_dna_bams = self.load_df_bams(is_dna=True)
        df_rna_bams = self.load_df_bams(is_dna=False)
        df_bams_tumor = df_dna_bams[df_dna_bams.is_tumor].merge(
            df_rna_bams, on="patient_id", how="outer", suffixes=["_dna", "_rna"])[[
                "patient_id", "bam_id_dna", "bam_id_rna",
                "sample_id_rna", "sample_id_dna"]]
        df_bams_tumor.rename(columns={"bam_id_dna": "bam_id_dna_tumor",
                                      "bam_id_rna": "bam_id_rna_tumor",
                                      "sample_id_dna": "sample_id_dna_tumor",
                                      "sample_id_rna": "sample_id_rna_tumor"},
                             inplace=True)
        df_bams_normal = df_dna_bams[~df_dna_bams.is_tumor].reset_index(drop=True)[["patient_id", "bam_id", "sample_id"]]
        df_bams_normal.rename(columns={"bam_id": "bam_id_dna_normal", "sample_id": "sample_id_dna_normal"}, inplace=True)
        df_bams = df_bams_normal.merge(df_bams_tumor, on="patient_id")

        df_clinical = self.load_df_clinical()
        pfs_col = self.get_pfs_col()

        df = df_clinical.merge(df_bams, on="patient_id", how="outer")

        df_hla = self.load_df_hla()
        len_df = len(df)
        df = df.merge(df_hla, on="patient_id", how="outer")
        assert len(df) == len_df, "Length changed when merging with HLA"

        patients = []
        for i, row in df.iterrows():
            patient_id = row["patient_id"]
            snv_vcf_paths = []
            indel_vcf_paths = []
            normal_sample = None
            tumor_sample = None
            if type(row["bam_id_dna_normal"]) == str:
                bam_id_dna_normal = row["bam_id_dna_normal"]
                bam_id_dna_tumor = row["bam_id_dna_tumor"]
                bam_id_rna_tumor = row["bam_id_rna_tumor"]
                bam_suffix = ".bam" if self.without_bqsr else ".bam-indelrealigned-default-default-bqsr-default-default.bam"
                bam_path_dna_normal = path.join(self.bam_dna_data_dir,
                                                bam_id_dna_normal) + bam_suffix
                bam_path_dna_tumor = path.join(self.bam_dna_data_dir,
                                               bam_id_dna_tumor) + bam_suffix
                bam_path_rna_tumor = (path.join(self.bam_rna_data_dir,
                                               bam_id_rna_tumor) +
                                      "-b2fq-PE-star-alignedAligned.sortedByCoord.out.bam")
                kallisto_dirname = "{}-kallisto".format(bam_id_rna_tumor)
                kallisto_path = path.join(self.kallisto_data_dir,
                                          kallisto_dirname,
                                          "abundance.tsv")
                cufflinks_dirname = "%s-b2fq-PE-star-alignedAligned.sortedByCoord.out.bam-cufflinks_output" % bam_id_rna_tumor
                cufflinks_path = path.join(self.cufflinks_data_dir,
                                           cufflinks_dirname,
                                           "genes.fpkm_tracking")
                normal_sample = Sample(
                    is_tumor=False,
                    bam_path_dna=bam_path_dna_normal)
                tumor_sample = Sample(
                    is_tumor=True,
                    bam_path_dna=bam_path_dna_tumor,
                    bam_path_rna=bam_path_rna_tumor,
                    kallisto_path=kallisto_path,
                    cufflinks_path=cufflinks_path)
                if self.without_bqsr:
                    snv_vcf_paths=[
                        path.join(self.vcf_data_dir, mutect_without_bqsr_snv_file_format_func(
                            patient_id, bam_id_dna_normal, bam_id_dna_tumor)),
                        path.join(self.vcf_data_dir, strelka_without_bqsr_snv_file_format_func(
                            patient_id, bam_id_dna_normal, bam_id_dna_tumor))]
                    indel_vcf_paths = []
                else:
                    snv_vcf_paths=[
                        path.join(self.vcf_data_dir, mutect_snv_file_format_func(
                            patient_id, bam_id_dna_normal, bam_id_dna_tumor)),
                        path.join(self.vcf_data_dir, strelka_snv_file_format_func(
                            patient_id, bam_id_dna_normal, bam_id_dna_tumor))]
                    indel_vcf_paths=[
                        path.join(self.vcf_data_dir, strelka_indel_file_format_func(
                            patient_id, bam_id_dna_normal, bam_id_dna_tumor))]
            elif self.only_patients_with_bams:
                # Don't create a Patient for non-bam patients if only_patients_with_bams
                # is True.
                continue

            if patient_id not in self.exclude_patient_ids:
                patient = Patient(id=patient_id,
                                  benefit=row["is_benefit"],
                                  os=row["OS in days"],
                                  pfs=row[pfs_col], # Depends in RECIST choice
                                  deceased=row["is_deceased"],
                                  progressed=row["is_progressed"],
                                  progressed_or_deceased=row["is_progressed_or_deceased"],
                                  hla_alleles=row["hla_allele_list"],
                                  snv_vcf_paths=snv_vcf_paths,
                                  indel_vcf_paths=indel_vcf_paths,
                                  normal_sample=normal_sample,
                                  tumor_sample=tumor_sample,
                                  additional_data=row.to_dict())
                patients.append(patient)

        extra_df_loaders = []
        extra_df_loaders.append(DataFrameLoader(
            "tcr", self.load_df_tcr))
        extra_df_loaders.append(DataFrameLoader(
            "tcr_peripheral_a", self.load_df_tcr_peripheral_a))
        extra_df_loaders.append(DataFrameLoader(
            "tcr_tumor", self.load_df_tcr_tumor))
        extra_df_loaders.append(DataFrameLoader(
            "pdl1", self.load_df_pdl1))
        extra_df_loaders.append(DataFrameLoader(
            "pdl1_all", lambda: self.load_df_pdl1(all_locations=True)))
        extra_df_loaders.append(DataFrameLoader(
            "cibersort", self.load_df_cibersort))
        extra_df_loaders.append(DataFrameLoader(
            "tcr_raw", self.load_df_tcr_raw))
        extra_df_loaders.append(DataFrameLoader(
            "tcr_expansion_a_b", self.load_df_tcr_expansion_a_b))
        extra_df_loaders.append(DataFrameLoader(
            "tcr_expansion_a_c", self.load_df_tcr_expansion_a_c))
        extra_df_loaders.append(DataFrameLoader(
            "genentech_tcga_subtypes", self.load_df_genentech_tcga_subtypes))
        extra_df_loaders.append(DataFrameLoader(
            "learned_genentech_tcga_subtypes", self.load_df_learned_tcga_subtypes))
        extra_df_loaders.append(DataFrameLoader(
            "signatures", self.load_df_signatures))

        cohort = Cohort(patients,
                        kallisto_ensembl_version=75,
                        cache_dir=self.cache_data_dir,
                        extra_df_loaders=extra_df_loaders,
                        join_with=self.join_with,
                        join_how=self.join_how,
                        filter_fn=self.filter_fn,
                        normalized_per_mb=self.normalized_per_mb,
                        min_coverage_normal_depth=self.min_coverage_normal_depth,
                        min_coverage_tumor_depth=self.min_coverage_tumor_depth,
                        responder_pfs_equals_os=self.responder_pfs_equals_os,
                        polyphen_dump_path=self.polyphen_data_dump,
                        pageant_coverage_path=self.pageant_coverage_data_dir,
                        benefit_plot_name=self.benefit_plot_name,
                        **self.init_kwargs)

        # Add some extra variables
        cohort.benefit_os_plot_name = self.benefit_os_plot_name
        cohort.hazard_plot_name = self.hazard_plot_name
        cohort.hazard_os_plot_name = self.hazard_os_plot_name
        cohort.tcr_peripheral_a_clonality_plot_name = self.tcr_peripheral_a_clonality_plot_name
        cohort.tcr_peripheral_a_clonality_short_plot_name = self.tcr_peripheral_a_clonality_short_plot_name
        cohort.tcr_tumor_clonality_plot_name = self.tcr_tumor_clonality_plot_name
        cohort.til_fraction_plot_name = self.til_fraction_plot_name

        patient_ids_without_bams = {"0979", "7592", "8214"}
        if self.only_patients_with_bams:
            num_expected_patients = 29 - len(self.exclude_patient_ids.union(patient_ids_without_bams))
        else:
            num_expected_patients = 29 - len(self.exclude_patient_ids)
        assert len(cohort) == num_expected_patients, (
            "Expected %d patients but have %d" % (num_expected_patients,
                                                  len(cohort)))
        return cohort

    def get_pfs_col(self):
        return "PFS (RECIST 1.1) in days"

    def load_df_clinical(self):
        # Add clinical data
        df_clinical = pd.read_csv(path.join(self.repo_data_dir, self.clinical_csv_filename), encoding="latin1")
        for col in df_clinical.columns:
            df_clinical.rename(columns={col: col.strip()}, inplace=True)
        df_clinical.rename(columns={"ID": "patient_id"}, inplace=True)
        df_clinical.patient_id = df_clinical.patient_id.apply(patient_id_to_str)

        os_col = "OS in days"
        pfs_col = self.get_pfs_col()

        # From NEJM, for reference: "Clinical benefit was defined as stable disease or better for
        # greater than 24 weeks after initiation  therapy"
        # Here we define it as 182 days.
        df_clinical["is_benefit"] = df_clinical[pfs_col].map(
            lambda pfs: pfs > np.floor(365 / 2.0))
        df_clinical["is_benefit_os"] = df_clinical[os_col].map(
            lambda os: os > 365)
        df_clinical["is_late_deceased"] = df_clinical[os_col].map(
            lambda os: os > np.floor(365 / 4.0))
        df_clinical["is_deceased"] = df_clinical["Alive Status"] != "Y"

        # Cannot use pfs < os as the definition of "is_progressed"; instead,
        # use the ongoing responder column.
        df_clinical["is_progressed"] = df_clinical["Ongoing Responder RECIST 1.1"].map(lambda y_or_n: y_or_n == "N")
        df_clinical["is_progressed_or_deceased"] = df_clinical.is_progressed | df_clinical.is_deceased

        # Sanity checking
        deceased_ids = set(df_clinical[df_clinical.is_deceased].patient_id.unique())
        progressed_or_deceased_ids = set(
            df_clinical[df_clinical.is_progressed_or_deceased].patient_id.unique())
        # Deceased is a subset of progressed/deceased
        assert progressed_or_deceased_ids.intersection(deceased_ids) == deceased_ids
        # No benefit is a subset of progressed; all no benefitters progressed
        # (and some benefitters progressed later on)
        progressed_ids = set(
            df_clinical[df_clinical.is_progressed].patient_id.unique())
        no_benefit_ids = set(df_clinical[~df_clinical.is_benefit].patient_id.unique())
        assert no_benefit_ids.intersection(progressed_ids) == no_benefit_ids

        return df_clinical

    def load_df_bams(self, is_dna=True):
        # Use the sequencing manifest to link sample IDs with BAMs
        df_sequencing = pd.read_csv(path.join(self.repo_data_dir, "sequencing_manifest.csv"))
        df_sequencing.rename(columns={"Individual Name": "patient_id", "Sample Name": "sample_id"}, inplace=True)
        df_sequencing.patient_id = df_sequencing.patient_id.apply(patient_id_to_str)

        # Sample 5122 has a comma in the MGI DNA name that shouldn't be there for the next join
        df_sequencing["MGI DNA Name"] = df_sequencing["MGI DNA Name"].apply(lambda name: "".join(name.split(",")))

        manifest_filename = "bam_manifest.csv" if is_dna else "2850417_Neoantigen_RNA_bams.csv"
        df_bams = pd.read_csv(path.join(self.repo_data_dir, manifest_filename), delimiter="\t", header=None)

        # Use the BAM manifest to get BAM IDs
        bam_id_col = 4 if is_dna else 1
        df_bams.rename(columns={0: "MGI DNA Name", bam_id_col: "bam_id"}, inplace=True)
        df_bam_map = df_bams.merge(df_sequencing, on="MGI DNA Name", how="inner")
        df_bam_map.bam_id = df_bam_map.bam_id.apply(lambda id: path.basename(id).split(".bam")[0])

        # Tumor vs. normal BAM
        include_columns = ["patient_id", "sample_id", "bam_id"]
        if is_dna:
            df_bam_map["is_tumor"] = df_bam_map["Common Name"].apply(lambda name: "normal" not in name)
            include_columns.append("is_tumor")
        df_bam_map = df_bam_map[include_columns]

        # Remove excluded samples
        # Sometimes we might try to exclude a sample that wasn't in the BAM list at all; use patient_ids_actually_excluded
        # to properly calculate the number of expected resultant rows.
        patient_ids_actually_excluded = self.exclude_patient_ids.intersection(set(df_bam_map.patient_id.unique()))
        if len(patient_ids_actually_excluded) > 0:
            df_bam_map = df_bam_map[~df_bam_map.patient_id.isin(patient_ids_actually_excluded)]
        expected_num = 26 - len(patient_ids_actually_excluded)
        if is_dna:
            expected_num *= 2
        assert len(df_bam_map) == expected_num, "Expected %d samples but only found %d" % (expected_num, len(df_bam_map))
        return df_bam_map

    def load_df_hla(self):
        # Grab the HLA alleles as a list of lists, in sample ID order
        df_hla = pd.read_csv(path.join(self.repo_data_dir, "hla", "optitype.tsv"), delimiter="\t")
        df_hla.rename(columns={"id": "bam_id"}, inplace=True)

        df_bams = self.load_df_bams(is_dna=True)
        df_hla = df_hla.merge(df_bams, on="bam_id")
        df_hla = df_hla[~df_hla.is_tumor]
        num_expected_allele_lists = len(df_bams) / 2
        assert len(df_hla) == num_expected_allele_lists, "Expected %d allele lists but only found %d" % (num_expected_allele_lists, len(df_hla))
        df_hla.sort_values(by="patient_id", inplace=True)
        def get_allele_list(row):
            return [row["A1"], row["A2"], row["B1"], row["B2"], row["C1"], row["C2"]]
        df_hla["hla_allele_list"] = df_hla.apply(get_allele_list, axis=1)
        return df_hla[["patient_id", "hla_allele_list"]]

    def load_df_tcr(self):
        # TCR-Seq clonality at peripheral timepoint A
        df_tcr = pd.read_csv(path.join(self.repo_data_dir, "tcr_master.csv"),
                             dtype={"Subject ID": object})
        df_tcr.rename(columns={"Subject ID": "patient_id"}, inplace=True)
        df_tcr.patient_id = df_tcr.patient_id.apply(lambda s: s.split("-")[0])
        return df_tcr

    def load_df_tcr_raw(self):
        df_tcr = self.load_df_tcr()
        dfs_tcr_raw = []
        for i, row in df_tcr.iterrows():
            sample_id = row["Sample ID"]
            df_tcr_raw = pd.read_csv(path.join(
                self.tcr_data_dir, sample_id + ".tsv"),
                sep="\t", usecols=["amino_acid", "rearrangement", "frequency"])
            df_tcr_raw["patient_id"] = row["patient_id"]
            df_tcr_raw["time_point"] = row["Time Point"]
            df_tcr_raw["sample_type"] = row["Sample Type"]
            df_tcr_raw = df_tcr_raw[["patient_id", "time_point", "sample_type", "amino_acid", "rearrangement", "frequency"]]
            dfs_tcr_raw.append(df_tcr_raw)
        return pd.concat(dfs_tcr_raw)

    def load_df_tcr_peripheral_a(self):
        df_tcr = self.load_df_tcr()
        df_tcr_peripheral_a = df_tcr[(df_tcr["Time Point"] == "A") &
                                     (df_tcr["Sample Type"] == "PBMC")]
        return df_tcr_peripheral_a

    def load_df_tcr_tumor(self):
        df_tcr = self.load_df_tcr()
        df_tcr_tumor = df_tcr[df_tcr["Sample Type"] == "Tumor"]
        return df_tcr_tumor

    def load_df_pdl1(self, all_locations=False):
        # Get PD-L1 expression data
        df_pdl1 = pd.read_csv(path.join(self.repo_data_dir, "pdl1_ihc.csv"))
        df_pdl1_id_map = pd.read_csv(path.join(self.repo_data_dir, "genentech_to_msk_id_map.csv"))
        df_pdl1 = df_pdl1.merge(df_pdl1_id_map, left_on="Patient Enrolled ID", right_on="Genentech Pt ID")
        df_pdl1["Sample ID"] = df_pdl1["Sample ID"].apply(lambda s: s.replace("-", ""))
        df_pdl1.rename(columns={"Sample ID": "patient_id"}, inplace=True)
        # Fix an ID mistake: 2397 in PD-L1 data = 2937 elsewhere
        df_pdl1.patient_id = df_pdl1.patient_id.apply(lambda sample: sample if sample != "2397" else "2937")
        # Manually selected (by Alex S.) locations for PD-L1 IHC values when patients
        # have multiple samples
        patient_locations = {
            "1022": "BLADDER",
            "1023": "BLADDER",
            "1026": "BLADDER RADICAL CYSTECTOMY",
            "1184": "LYMPHNODES",
            "1202": "BLADDER URETUS BILATERAL TUBES",
            "1232": "BLADDER"
        }

        def filter_pdl1_rows(df, patient_locations):
            kept_rows = []
            for i, row in df.iterrows():
                genentech_id = str(row["Patient Enrolled ID"])
                if genentech_id in patient_locations:
                    location = row["Location"]
                    if location in patient_locations[genentech_id]:
                        kept_rows.append(row)
                else:
                    kept_rows.append(row)
            return pd.DataFrame.from_records(kept_rows, columns=df.columns)

        if all_locations:
            return df_pdl1
        return filter_pdl1_rows(df_pdl1, patient_locations)

    def load_df_cibersort(self):
        df_cibersort = pd.read_csv(
            path.join(self.repo_data_dir, "bladder-cibersort.txt"), sep="\t")
        df_bams = self.load_df_bams(is_dna=False)
        df_cibersort["bam_id"] = df_cibersort["SampleId"].str.split("/").str.get(1).str.rsplit("-", n=5).str.get(0)
        df_cibersort = df_cibersort.merge(df_bams, on="bam_id")
        return df_cibersort

    def load_df_tcr_expansion(self, time_point):
        df_tcr_expansion = pd.read_csv(path.join(self.repo_data_dir, "tcr_expansion.csv"))
        # Some rows have NaNs for # of expanded clones and # of expanded TIL clones; remove those rows.
        df_tcr_expansion = df_tcr_expansion.dropna(subset=["N Expanded Clones"])
        df_tcr_expansion["patient_id"] = df_tcr_expansion["Subject ID"].apply(lambda s: s.split("-")[0])
        df_tcr_expansion = df_tcr_expansion[df_tcr_expansion["Time Point"] == time_point]
        return df_tcr_expansion

    def load_df_tcr_expansion_a_b(self):
        return self.load_df_tcr_expansion("B")

    def load_df_tcr_expansion_a_c(self):
        return self.load_df_tcr_expansion("C")

    def load_df_genentech_tcga_subtypes(self):
        df_tcga_subtypes = pd.read_csv(path.join(self.repo_data_dir, "tcga_subtypes.csv"), header=1)
        df_tcga_subtype_ids = pd.read_csv(path.join(self.repo_data_dir, "tcga_subtypes_id_map.csv"), header=0)
        # Some IDs look like e.g. 00-40
        df_tcga_subtype_ids["patient_id"] = df_tcga_subtype_ids["Sample ID"].apply(
            lambda s: "".join([ch for ch in s if ch.isdigit()]))
        df_tcga_subtypes_merged = df_tcga_subtype_ids.merge(df_tcga_subtypes,
                                                    left_on="Genentech Pt ID",
                                                    right_on="Patient Enrolled ID")
        df_tcga_subtypes_merged = df_tcga_subtypes_merged[~df_tcga_subtypes_merged["TCGA Subtype"].isin(
            [".", "Unclassified"])]
        df_tcga_subtypes_merged["basal"] = df_tcga_subtypes_merged["TCGA Subtype"].apply(lambda x: False if x in ["I", "II"] else True)

        return df_tcga_subtypes_merged

    def load_df_learned_tcga_subtypes(self):
        df_tcga_subtypes = \
            pd.read_csv(path.join(self.repo_data_dir, "all_tcga_subtypes.csv"), dtype = {"patient_id" : str})
        df_tcga_subtypes["basal"] = df_tcga_subtypes["TCGA Subtype"].apply(lambda x: False if x in [1,2] else True)
        return df_tcga_subtypes

    def load_df_signatures(self):
        sigs_df = pd.read_csv(path.join(self.repo_data_dir, "rizvi_cohort_signatures.csv"),
                              dtype={"Sample": str})
        sigs_df.rename(columns={"Sample": "patient_id"}, inplace=True)
        return sigs_df

def patient_id_to_str(patient_id):
    return str(patient_id).rjust(4, "0")
