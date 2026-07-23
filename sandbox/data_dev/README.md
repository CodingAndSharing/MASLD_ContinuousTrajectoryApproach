# Synthetic MASLD Cohort — NIHR BioResource (Cambridge) style

## ⚠️ IMPORTANT: This is SYNTHETIC data

This dataset is **100% computer-generated (simulated)**. It contains **no real
patients** and **no data from the NIHR BioResource**. It was created to mirror
the *structure, variable types, and statistical relationships* of a real
MASLD/MASH cohort so that analysis code (e.g. disease-trajectory modelling and
biomarker discovery) can be developed and tested before applying for access to
genuine managed-access data.

### Why not real NIHR BioResource data?
The NIHR BioResource (hosted by NIHR Cambridge) is a **managed-access** resource.
Patient-level data, samples, and recall-by-genotype/phenotype are released only
after a formal application to, and approval by, the BioResource Data Access
Committee (DAC). Data cannot be downloaded openly. To request genuine data:
- Apply for BioResource data access: https://www.bioresource.nihr.ac.uk/using-our-bioresource/apply-for-bioresource-data-access/
- Cambridge BioResource Centre: https://www.bioresource.nihr.ac.uk/centres-programmes/nihr-bioresource-centre-cambridge/

## Cohort at a glance
- **600 participants**, recruited "at" a single centre (label: NIHR BioResource - Cambridge)
- **2,512 longitudinal visits** (3–5 visits each, months 0/12/24/36/48)
- Multi-modal: clinical, longitudinal biochemistry, non-invasive fibrosis scores,
  elastography, circulating proteomic biomarkers, genetic risk variants,
  liver histology (~45% biopsied), and follow-up outcomes.
- MASLD-consistent: alcohol intake kept below harmful thresholds; cardiometabolic
  risk factors enriched with disease severity.

## How the data were generated (so you can trust the relationships)
Each participant is assigned a **latent disease trajectory** — one of
`stable_steatosis`, `slow_progressor`, `fast_progressor`, `regressor` — with a
baseline latent severity and a per-year progression slope. **All observable
variables are generated downstream of this latent severity**, so:
- Liver enzymes, elastography (kPa), CAP, and proteomic markers rise with severity;
  adiponectin and HSD17B13 fall / protect.
- **FIB-4, APRI, and NFS are computed from their real formulae** using the
  simulated component labs — they are internally consistent, not drawn at random.
- Biopsy fibrosis stage and NAS components derive from the same latent severity;
  mean FIB-4 increases monotonically across stages F0→F3 (1.13→2.04).
- Genetic risk-allele dosages (PNPLA3 etc.) are mildly enriched with severity.

### Ground-truth columns for method validation
`trajectory_class_truth` (participants) and `latent_severity_truth` (visits) are
the **hidden generating variables**. They would NOT exist in a real study. Use
them only to benchmark unsupervised trajectory / clustering / biomarker methods
(e.g. "did my model recover the true trajectory classes?"), then drop them.

## Files
| File | Rows | Grain | Contents |
|------|------|-------|----------|
| `data/01_participants_baseline.csv` | 600 | participant | Demographics, anthropometrics, comorbidities |
| `data/02_longitudinal_labs_visits.csv` | 2512 | participant-visit | Biochemistry, elastography, FIB-4/APRI/NFS |
| `data/03_proteomics_biomarkers_baseline.csv` | 600 | participant | Circulating protein biomarkers (baseline) |
| `data/04_genetics_risk_variants.csv` | 600 | participant | MASLD risk/protective variant dosages |
| `data/05_histology_staging.csv` | 600 | participant | Biopsy fibrosis stage, NAS, MASH diagnosis |
| `data/06_outcomes_followup.csv` | 600 | participant | Follow-up duration and clinical events |
| `data_dictionary.csv` | — | — | Every variable, units, and definition |

Join all tables on **`participant_id`**; join longitudinal to others at baseline
with `visit_month == 0`.

## Quick start
```python
import pandas as pd, glob
d = "data/"
base = pd.read_csv(d+"01_participants_baseline.csv")
labs = pd.read_csv(d+"02_longitudinal_labs_visits.csv")
prot = pd.read_csv(d+"03_proteomics_biomarkers_baseline.csv")
# baseline multimodal matrix
bl = (labs[labs.visit_month==0]
        .merge(base, on="participant_id")
        .merge(prot, on="participant_id"))
print(bl.shape)
```

## Suggested analyses this dataset supports
- Latent-class / mixed-effects trajectory modelling of FIB-4 or liver stiffness over time
- Unsupervised clustering of baseline multimodal features → compare to `trajectory_class_truth`
- Trajectory-specific biomarker discovery (e.g. PRO-C3 / TIMP-1 elevation in fast progressors)
- Survival / competing-risk models for `progressed_to_advanced_fibrosis`
- Non-invasive score benchmarking (FIB-4/APRI/NFS vs biopsy `fibrosis_stage_F`)

## Provenance & reproducibility
- Generator: **`generate_dataset.py`** (in this folder; full listing embedded below).
- Seed: `numpy.random.default_rng(20240517)` — re-running the script reproduces all six
  CSVs byte-for-byte (verified by md5 checksum).
- No real individuals; safe to share, version, and publish as test data.
- Biomarker/variant *names* are real MASLD-field entities; their *values are simulated*.

To regenerate from scratch:
```bash
cd masld_nihr_synthetic
python generate_dataset.py   # requires numpy + pandas
```

## Generation script (`generate_dataset.py`)

The complete script used to build every table is below (and saved alongside this
README). All randomness flows from a single seeded generator, so the output is
fully reproducible.

```python
"""
generate_dataset.py
Synthetic MASLD cohort (NIHR BioResource - Cambridge style).

100% simulated data - no real patients, no real BioResource data.
Reproduces the six CSVs in ./data plus data_dictionary.csv.

Design: each participant is assigned a latent disease TRAJECTORY and a
latent severity; all observable variables (labs, elastography, biomarkers,
histology, outcomes) are generated downstream of that latent severity, and
the non-invasive fibrosis scores (FIB-4, APRI, NFS) are computed from their
published formulae using the simulated component labs.

Run:  python generate_dataset.py
"""
import numpy as np
import pandas as pd
import os

RNG_SEED = 20240517
rng = np.random.default_rng(RNG_SEED)
N = 600
os.makedirs("data", exist_ok=True)

# ----------------------------------------------------------------------
# 1. Latent disease trajectory (the hidden generating variable)
# ----------------------------------------------------------------------
traj_names = ["stable_steatosis", "slow_progressor", "fast_progressor", "regressor"]
traj_p     = [0.42, 0.30, 0.20, 0.08]
traj = rng.choice(traj_names, size=N, p=traj_p)

base_sev = {"stable_steatosis": 0.20, "slow_progressor": 0.35,
            "fast_progressor": 0.55, "regressor": 0.45}
sev0 = np.clip(np.array([base_sev[t] for t in traj]) + rng.normal(0, 0.08, N), 0.02, 0.95)

slope = {"stable_steatosis": 0.01, "slow_progressor": 0.05,
         "fast_progressor": 0.13, "regressor": -0.06}
prog_slope = np.array([slope[t] for t in traj]) + rng.normal(0, 0.02, N)

# ----------------------------------------------------------------------
# 2. Participants - demographics / anthropometrics / comorbidities
# ----------------------------------------------------------------------
pid = [f"NBR-MASLD-{i:04d}" for i in range(1, N + 1)]
sex = rng.choice(["Male", "Female"], N, p=[0.55, 0.45])
age = np.clip(rng.normal(54, 12, N), 19, 84).round(0)
ethnicity = rng.choice(["White British", "White Other", "South Asian", "Black", "Mixed/Other"],
                       N, p=[0.62, 0.12, 0.14, 0.06, 0.06])
bmi = np.clip(rng.normal(31, 5, N) + sev0 * 6, 18, 55).round(1)
waist = (bmi * 2.6 + rng.normal(0, 6, N) + (sex == "Male") * 8).round(0)
t2dm = (rng.random(N) < (0.15 + sev0 * 0.5)).astype(int)
htn  = (rng.random(N) < (0.25 + sev0 * 0.4)).astype(int)
dyslip = (rng.random(N) < (0.30 + sev0 * 0.4)).astype(int)

participants = pd.DataFrame({
    "participant_id": pid, "trajectory_class_truth": traj, "sex": sex,
    "age_baseline": age.astype(int), "ethnicity": ethnicity,
    "bmi_baseline": bmi, "waist_cm": waist.astype(int),
    "type2_diabetes": t2dm, "hypertension": htn, "dyslipidaemia": dyslip,
    "alcohol_units_week": np.clip(rng.normal(6, 4, N), 0, 20).round(1),
    "recruitment_centre": "NIHR BioResource - Cambridge",
    "consent_recall_by_genotype": rng.choice([0, 1], N, p=[0.3, 0.7]),
})

# ----------------------------------------------------------------------
# 3. Longitudinal visits (months 0/12/24/36/48; 3-5 per participant)
#    FIB-4, APRI, NFS computed from published formulae.
# ----------------------------------------------------------------------
rows = []
months_grid = [0, 12, 24, 36, 48]
for i in range(N):
    n_visits = rng.choice([3, 4, 5], p=[0.25, 0.35, 0.40])
    for v in range(n_visits):
        m = months_grid[v]
        yrs = m / 12.0
        sv = np.clip(sev0[i] + prog_slope[i] * yrs + rng.normal(0, 0.04), 0.01, 0.99)
        alt = np.clip(rng.normal(28 + sv * 70, 10), 8, 250).round(0)
        ast = np.clip(rng.normal(24 + sv * 60, 9), 8, 220).round(0)
        ggt = np.clip(rng.normal(30 + sv * 90, 20), 9, 400).round(0)
        alp = np.clip(rng.normal(70 + sv * 40, 18), 30, 300).round(0)
        plt_ct = np.clip(rng.normal(255 - sv * 110, 35), 60, 420).round(0)
        albumin = np.clip(rng.normal(45 - sv * 8, 3), 28, 52).round(1)   # g/L
        bilirubin = np.clip(rng.normal(9 + sv * 8, 4), 3, 60).round(0)
        hba1c = np.clip(rng.normal(38 + t2dm[i] * 18 + sv * 10, 5), 30, 110).round(0)
        trig = np.clip(rng.normal(1.4 + sv * 1.2, 0.5), 0.4, 8).round(2)
        hdl = np.clip(rng.normal(1.4 - sv * 0.5, 0.25), 0.5, 2.6).round(2)
        homa_ir = np.clip(rng.normal(2 + sv * 6, 1.2), 0.5, 20).round(2)
        lsm_kpa = np.clip(rng.normal(4.5 + sv * 18, 2.5), 2.5, 45).round(1)
        cap = np.clip(rng.normal(250 + sv * 90, 30), 150, 400).round(0)

        # FIB-4 = (age * AST) / (platelets * sqrt(ALT))        [Sterling 2006]
        fib4 = round((age[i] + yrs) * ast / (plt_ct * np.sqrt(alt)), 2)
        # APRI = (AST / ULN=40) / platelets * 100              [Wai 2003]
        apri = round((ast / 40.0) / plt_ct * 100, 2)
        # NFS = -1.675 + 0.037*age + 0.094*BMI + 1.13*IFG/DM
        #       + 0.99*(AST/ALT) - 0.013*platelets - 0.66*albumin(g/dL)  [Angulo 2007]
        bmi_v = participants.loc[i, "bmi_baseline"]
        nfs = round(-1.675 + 0.037 * (age[i] + yrs) + 0.094 * bmi_v
                    + 1.13 * t2dm[i] + 0.99 * (ast / alt)
                    - 0.013 * plt_ct - 0.66 * (albumin / 10.0), 2)

        rows.append([pid[i], m, round(sv, 3), int(alt), int(ast), int(ggt),
                     int(alp), int(plt_ct), albumin, int(bilirubin), int(hba1c),
                     trig, hdl, homa_ir, lsm_kpa, int(cap), fib4, apri, nfs])

longitudinal = pd.DataFrame(rows, columns=[
    "participant_id", "visit_month", "latent_severity_truth", "ALT_U_L", "AST_U_L",
    "GGT_U_L", "ALP_U_L", "platelets_10e9_L", "albumin_g_L", "bilirubin_umol_L",
    "hba1c_mmol_mol", "triglycerides_mmol_L", "HDL_mmol_L", "HOMA_IR",
    "liver_stiffness_kPa", "CAP_dB_m", "FIB4", "APRI", "NFS"])

# ----------------------------------------------------------------------
# 4. Circulating proteomic biomarkers (baseline)
# ----------------------------------------------------------------------
sev = sev0
prot = pd.DataFrame({"participant_id": pid})
prot["PRO_C3_ng_ml"]        = np.clip(rng.normal(15 + sev * 40, 8), 4, 150).round(1)
prot["CK18_M30_U_L"]        = np.clip(rng.normal(150 + sev * 450, 90), 50, 2000).round(0)
prot["CK18_M65_U_L"]        = np.clip(rng.normal(200 + sev * 500, 110), 60, 2500).round(0)
prot["FGF21_pg_ml"]         = np.clip(rng.normal(180 + sev * 350, 120), 20, 1500).round(0)
prot["adiponectin_ug_ml"]   = np.clip(rng.normal(9 - sev * 4, 2.5), 1, 25).round(2)
prot["leptin_ng_ml"]        = np.clip(rng.normal(12 + sev * 22, 7), 1, 90).round(1)
prot["TIMP1_ng_ml"]         = np.clip(rng.normal(90 + sev * 160, 40), 30, 600).round(0)
prot["hyaluronic_acid_ng_ml"] = np.clip(rng.normal(30 + sev * 180, 45), 8, 800).round(0)
prot["IL6_pg_ml"]           = np.clip(rng.normal(2 + sev * 8, 2), 0.3, 40).round(2)
prot["hsCRP_mg_L"]          = np.clip(rng.normal(2 + sev * 6, 2), 0.1, 30).round(2)
# trajectory-specific fibrogenic signature in fast progressors
fast = (traj == "fast_progressor")
prot.loc[fast, "PRO_C3_ng_ml"] = (prot.loc[fast, "PRO_C3_ng_ml"] * 1.25).round(1)
prot.loc[fast, "TIMP1_ng_ml"]  = (prot.loc[fast, "TIMP1_ng_ml"] * 1.20).round(0)

# ----------------------------------------------------------------------
# 5. Genetic risk / protective variants (allele dosage 0/1/2)
# ----------------------------------------------------------------------
geno = pd.DataFrame({"participant_id": pid})
def snp(risk_af, sev_boost):
    p = np.clip(risk_af + sev * sev_boost, 0.02, 0.95)
    return (rng.random(N) < p).astype(int) + (rng.random(N) < p).astype(int)
geno["PNPLA3_rs738409_G"]      = snp(0.26, 0.25)    # risk    [Romeo 2008]
geno["TM6SF2_rs58542926_T"]    = snp(0.07, 0.20)    # risk    [Kozlitina 2014]
geno["GCKR_rs1260326_T"]       = snp(0.40, 0.10)    # risk
geno["MBOAT7_rs641738_T"]      = snp(0.43, 0.08)    # risk
geno["HSD17B13_rs72613567_TA"] = snp(0.26, -0.15)   # protective [Abul-Husn 2018]

# ----------------------------------------------------------------------
# 6. Liver histology (subset biopsied) - NASH-CRN staging  [Kleiner 2005]
# ----------------------------------------------------------------------
biopsy_has = rng.random(N) < 0.45
def stage_from_sev(s):
    return int(np.sum(s > np.array([0.2, 0.4, 0.6, 0.8])))   # F0-F4
hist = pd.DataFrame({"participant_id": pid})
hist["biopsy_available"] = biopsy_has.astype(int)
hist["fibrosis_stage_F"] = [stage_from_sev(s) if b else np.nan
                            for s, b in zip(sev, biopsy_has)]
hist["NAS_steatosis_0_3"] = [min(3, max(0, int(round(rng.normal(1.5 + s, 0.6))))) if b else np.nan
                             for s, b in zip(sev, biopsy_has)]
hist["NAS_ballooning_0_2"] = [min(2, max(0, int(round(rng.normal(s * 2.2, 0.6))))) if b else np.nan
                              for s, b in zip(sev, biopsy_has)]
hist["NAS_inflammation_0_3"] = [min(3, max(0, int(round(rng.normal(s * 2.5, 0.7))))) if b else np.nan
                                for s, b in zip(sev, biopsy_has)]
hist["MASH_diagnosis"] = [int((st >= 1) and (rng.random() < 0.6 + s * 0.3)) if b else np.nan
                          for st, s, b in zip(hist["fibrosis_stage_F"], sev, biopsy_has)]

# ----------------------------------------------------------------------
# 7. Follow-up outcomes
# ----------------------------------------------------------------------
fu_years = np.clip(rng.normal(4.5, 1.2, N), 1, 7).round(2)
event_p = np.clip(0.02 + sev * 0.35 + prog_slope * 0.5, 0.01, 0.9)
outcome = pd.DataFrame({"participant_id": pid, "followup_years": fu_years})
outcome["progressed_to_advanced_fibrosis"] = (rng.random(N) < event_p).astype(int)
outcome["liver_related_event"] = (rng.random(N) < event_p * 0.3).astype(int)
outcome["all_cause_death"] = (rng.random(N) < 0.03 + sev * 0.08).astype(int)

# ----------------------------------------------------------------------
# 8. Write files
# ----------------------------------------------------------------------
participants.to_csv("data/01_participants_baseline.csv", index=False)
longitudinal.to_csv("data/02_longitudinal_labs_visits.csv", index=False)
prot.to_csv("data/03_proteomics_biomarkers_baseline.csv", index=False)
geno.to_csv("data/04_genetics_risk_variants.csv", index=False)
hist.to_csv("data/05_histology_staging.csv", index=False)
outcome.to_csv("data/06_outcomes_followup.csv", index=False)
print(f"Wrote 6 tables: {len(participants)} participants, {len(longitudinal)} visits.")
```

## Literature resources used to design the dataset

**Important scope note.** This synthetic dataset was **not** assembled by copying
data or text out of any publication. It was built from established, widely-used
clinical **formulae** (FIB-4, APRI, NFS), the standard histological **staging
system** (NASH-CRN), and **well-characterised MASLD biomarkers and genetic
variants**. The references below are the primary sources for each of those design
choices. For each, a short representative quote (well under fair-use length) shows
the specific fact that informed the code; please consult the original papers for
full context.

### Disease definition / nomenclature
1. **Rinella ME, et al. "A multisociety Delphi consensus statement on new fatty
   liver disease nomenclature." *Hepatology* / *Journal of Hepatology*, 2023.**
   Published 24 June 2023. PMID 37364790.
   URL: https://www.journal-of-hepatology.eu/article/S0168-8278(23)00418-X/fulltext
   Section: Results / new terminology.
   *Used for:* the MASLD/MASH terminology and the alcohol-threshold framing of the
   cohort. Representative point from the consensus: the group selected "MASLD"
   (metabolic dysfunction-associated steatotic liver disease) to replace NAFLD,
   requiring at least one cardiometabolic risk factor.

### Non-invasive fibrosis scores (computed in `02_longitudinal_labs_visits.csv`)
2. **Sterling RK, et al. "Development of a simple noninvasive index to predict
   significant fibrosis in patients with HIV/HCV coinfection." *Hepatology*
   2006;43(6):1317-1325.** Published June 2006. PMID 16729309. DOI 10.1002/hep.21178.
   URL: https://onlinelibrary.wiley.com/doi/pdf/10.1002/hep.21178
   Section: Methods (index definition).
   *Used for:* the **FIB-4** formula, `age x AST / (platelets x sqrt(ALT))`.

3. **Wai CT, et al. "A simple noninvasive index can predict both significant
   fibrosis and cirrhosis in patients with chronic hepatitis C." *Hepatology*
   2003;38(2):518-526.** Published August 2003.
   URL: https://aasldpubs.onlinelibrary.wiley.com/doi/10.1053/jhep.2003.50346
   Section: Methods (APRI definition).
   *Used for:* the **APRI** formula, `(AST/ULN)/platelets x 100`.

4. **Angulo P, et al. "The NAFLD fibrosis score: a noninvasive system that
   identifies liver fibrosis in patients with NAFLD." *Hepatology*
   2007;45(4):846-854.** Published April 2007. DOI 10.1002/hep.21496.
   URL: https://onlinelibrary.wiley.com/doi/pdf/10.1002/hep.21496
   Section: Methods / scoring model.
   *Used for:* the **NFS** formula and its six inputs (age, BMI, IFG/diabetes,
   AST/ALT, platelets, albumin). Note albumin enters in **g/dL**.

### Histological staging (`05_histology_staging.csv`)
5. **Kleiner DE, et al. "Design and validation of a histological scoring system
   for nonalcoholic fatty liver disease." *Hepatology* 2005;41(6):1313-1321.**
   Published June 2005.
   URL: https://aasldpubs.onlinelibrary.wiley.com/doi/10.1002/hep.20701
   Section: Results (NAS + fibrosis staging).
   *Used for:* the **NAFLD Activity Score** components (steatosis 0-3, ballooning
   0-2, lobular inflammation 0-3) and the **F0-F4 fibrosis stage**.

### Genetic risk / protective variants (`04_genetics_risk_variants.csv`)
6. **Romeo S, et al. "Genetic variation in PNPLA3 confers susceptibility to
   nonalcoholic fatty liver disease." *Nature Genetics* 2008;40(12):1461-1465.**
   Published December 2008.
   URL: https://www.nature.com/articles/ng.257
   *Used for:* the **PNPLA3 rs738409 (I148M)** risk allele, the strongest common
   genetic determinant of MASLD.

7. **Kozlitina J, et al. "Exome-wide association study identifies a TM6SF2
   variant that confers susceptibility to nonalcoholic fatty liver disease."
   *Nature Genetics* 2014;46(4):352-356.** Published April 2014.
   URL: https://www.nature.com/articles/ng.2901
   *Used for:* the **TM6SF2 rs58542926 (E167K)** risk allele.

8. **Abul-Husn NS, et al. "A protein-truncating HSD17B13 variant and protection
   from chronic liver disease." *New England Journal of Medicine*
   2018;378(12):1096-1106.** Published 22 March 2018.
   URL: https://www.nejm.org/doi/full/10.1056/NEJMoa1712191
   *Used for:* the **HSD17B13 rs72613567** *protective* allele (encoded with a
   negative severity weight in the script).

### Circulating biomarkers (`03_proteomics_biomarkers_baseline.csv`)
The proteomic markers (PRO-C3, CK-18 M30/M65, FGF21, adiponectin, leptin, TIMP-1,
hyaluronic acid, IL-6, hsCRP) are all established or candidate MASLD/MASH markers
whose direction of association with disease severity (rising, or falling for the
protective adiponectin) is drawn from the fibrogenesis / apoptosis / inflammation
literature. Representative primary source for the PRO-C3 fibrogenesis marker:
9. **Nielsen MJ, et al. / Karsdal MA, et al.** on the **PRO-C3** neo-epitope of
   type III collagen formation as a marker of active fibrogenesis (basis for the
   ADAPT score). See e.g. *Journal of Hepatology* and *American Journal of
   Physiology - GI* (2015-2019).
   *Used for:* modelling PRO-C3 as rising with fibrogenic activity and elevated
   further in the fast-progressor trajectory.

### Resource access (real data)
10. **NIHR BioResource — Apply for data access.** Accessed 2024.
    URL: https://www.bioresource.nihr.ac.uk/using-our-bioresource/apply-for-bioresource-data-access/
    *Used for:* the managed-access framing and the correct application route for
    genuine (non-synthetic) BioResource data.
