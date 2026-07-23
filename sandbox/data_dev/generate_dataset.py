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
