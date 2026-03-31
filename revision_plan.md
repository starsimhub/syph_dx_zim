# Revision Plan: Syphilis Dx Zimbabwe

**Date:** 2026-03-31
**Branch:** `submission-mar31`
**Source:** Co-author feedback call (Lori Newman, Darcy Rao, Remco Peters), March 23, 2026
**Manuscript:** `syph_dx_zim_v2.md` (working file — v2 has partial edits from Lori/Darcy already applied)

---

## Status legend

- ✅ Done
- 🔄 In progress
- ⏳ Pending / not yet started
- ❌ Decided: will not do (this revision)
- 🔁 Requires recalibration — batch with Group 2

---

## Group 1: Text-only changes

| # | Change | Location | Status | Notes | RS Notes |
|---|--------|----------|--------|-------|----------|
| 1.1 | Remove newborn POC (CS) scenario | Abstract, Methods Table 2, Results CS paragraph, Fig 4 caption | ⏳ | Decided by Lori & Darcy. Remove CS row from Table 2, remove CS results paragraph, update Fig 4 caption and panel C. | OK TO PROCEED. WILL ALSO REQUIRE UPDATES TO FIGURES 2–4. |
| 1.2 | Standardize "pathways" → "screening and diagnostic strategies" | Throughout | ⏳ | Check: Table 1 title, section headers, results text, abstract | OK TO PROCEED |
| 1.3 | Reframe Table 1 columns | Methods Table 1 | ⏳ | Rename sensitivity columns: "Sensitivity (active infection)" and "False-positive rate (past/no active infection)" | OK TO PROCEED |
| 1.4 | Remove AMR mentions | Intro para 2; Limitations | ⏳ | Intro: remove "possible contribution to AMR concerns". Limitations: remove "nor did we model potential AMR to benzathine penicillin". | OK TO PROCEED |
| 1.5 | Strengthen benzathine penicillin supply angle | Policy implications | ⏳ | Already present — move to lead the paragraph, make it more prominent | OK TO PROCEED |
| 1.6 | Add HIV context (1–2 sentences) | Introduction | ⏳ | Brief mention of Zimbabwe's generalized HIV epidemic and syphilis-HIV coinfection as setting context | OK TO PROCEED |
| 1.7 | Investigate + fix Fig 3 vs Fig 4 number discrepancy | Plot scripts + Results text | ⏳ | Both claim treatments/year but numbers differ substantially (~15K GUD in Fig 3 vs ~140K total in Fig 4). Needs debugging of plot scripts before text can be updated. | NOT JUST A TIMING ISSUE — FIG 4 CLAIMS TREATMENTS PER YEAR. NEEDS DEBUGGING. |
| 1.8 | Update "confirmatory test" terminology | Throughout | ⏳ | Use "screening followed by confirmatory serological testing". Keep "confirmatory" where unambiguous; be approach-specific where not. Per Peters: calling it confirmatory is fine; don't specify T/NT. | PROCEED. PETERS: "CONFIRMATORY" IS PROBABLY FINE; PREFER "SCREENING FOLLOWED BY CONFIRMATORY SEROLOGICAL TESTING". |
| 1.9 | Update WHO department name for Peters affiliation | Author affiliations | ✅ | Done in v2 | |
| 1.10 | Fix author middle initials | Author list | ✅ | Alina M. Muellenmeister, Romesh G. Abeysuriya — done in v2 | |
| 1.11 | Refocus primary framing on overtreatment reduction | Abstract, Conclusions, throughout | ⏳ | Lori: emphasize combined strategy + ANC bottleneck. Frame syndromic management assumptions as caveats not conclusions ("assuming syndromic management is implemented as per guidelines..."). | DECISION: LORI — REFOCUS ON OVERTREATMENT REDUCTION VIA POC DIAGNOSTICS. SYNDROMIC MANAGEMENT SHOULD BE A CAVEAT NOT A CONCLUSION. |
| 1.12 | Include undertreatment section | Results/Discussion | ❌ | Not including in this paper. Possible follow-up. | DECISION: NOT GOING TO INCLUDE. CONSIDER BUDGET IMPACT ANALYSIS INSTEAD (see 1.12b). |
| 1.12b | Budget impact / cost-thresholding analysis | New discussion section or appendix | ⏳ | Not a full CEA. Threshold Q: given cost of GUD consult + penicillin, how much can a POC test cost and still be cost-saving? Scope = drug + test costs only, not health outcomes (BIA not CEA). | POSSIBLE ADDITION — BUDGET IMPACT NOT CEA. SCOPE AND DECISION STILL NEEDED. |
| 1.13 | Revise chancre visibility (men) and care-seeking (women) | Methods text — update after recalibration (2.2) | 🔁 | Men: chancre visibility 50% → 80% (heterosexual anatomy context). Women: care-seeking rate +50% relative to current calibrated values. Update methods text once recalibration done. | REQUIRES RECALIBRATION — SEE 2.2. |
| 1.14 | GUD sensitivity 80% — keep + explicit justification | Methods | ⏳ | Add a sentence making clear this is an assumption with brief justification. | KEEP; NOTE AS ASSUMPTION. |
| 1.15 | Remove syphilis transmission results | Results (Fig 1E), Discussion | ❌ | Decided: cut from this paper. Remove Fig 1E panel or simplify Fig 1 to 4 panels. | NOT GOING TO INCLUDE IN THIS PAPER. |
| 1.16 | Newborn pathway assumptions — clarify | Methods | ✅ | Resolved by removal of CS diagnostic scenario (1.1). | REMOVING THE NEWBORN CS DX ADDRESSES THIS. |

---

**COMMENTS FROM REMCO PETERS (on 1.8):**
Calling the second test "confirmatory" is probably fine. Prefer "screening followed by confirmatory serological testing" rather than specifying T or NT, since the confirmatory test could be either depending on the algorithm. For novel active-infection tests where the second test assigns infection status rather than improving specificity, be approach-specific.

**COMMENT FROM LORI (on 1.11):**
Refocus paper very clearly on overtreatment reduction via POC diagnostics (GUD + asymptomatic screening in high-risk groups + ANC). Be careful about "syndromic management works reasonably well" — this is an assumption, not a conclusion. Frame as: "assuming syndromic management is implemented as per guidelines..."

---

## Group 2: Changes requiring recalibration

**Do these as a single batch.** After recalibration: re-run top-200 msim and scenario analysis, then update methods text (1.13) and results.

| # | Change | Code to update | Status | Notes |
|---|--------|---------------|--------|-------|
| 2.1 | Update early/late latent cutoff to WHO 24 months | `diseases.py` — latent stage transition ~12 → 24 months | ⏳ | Will shift MTCT stage distribution (late latent currently 62% of adverse outcomes). Check sensitivity before committing; may affect overtreatment rates in screening pathways. |
| 2.2 | Increase chancre visibility (men) and care-seeking (women) | `interventions.py` (or `diseases.py`) — `p_symp` / care-seeking parameters | ⏳ | Men: chancre visibility 50% → 80%. Women: CS rate × 1.5 relative to current calibrated values. |
| 2.3 | Any follow-on parameter changes | `interventions.py`, `diseases.py` | ⏳ | Placeholder — check after all decisions resolved. |

---

## Group 3: Decisions resolved

| Topic | Decision | Notes |
|-------|----------|-------|
| Include undertreatment | ❌ Not in this paper | Possible follow-up; budget impact thresholding may be added (1.12b) |
| Newborn CS diagnostic scenario | ❌ Remove from paper | Lori & Darcy; resolved by 1.1 |
| "Confirmatory test" renaming | ✅ Keep "confirmatory"; use "screening followed by confirmatory serological testing" | Peters |
| Care-seeking / chancre visibility | 🔁 Revise with recalibration | Men 50→80% visibility; women CS ×1.5 |
| GUD sensitivity 80% | ✅ Keep; note as assumption | 1.14 |
| Syphilis transmission results | ❌ Remove from paper | 1.15 |

---

## Notes on v2 edits already made

The file `syph_dx_zim_v2.md` contains partial edits from Lori/Darcy:
- Title: "novel syphilis diagnostics" → "novel tests for active syphilis"
- Abstract: minor phrasing, "POC" abbreviation introduced
- Intro: GUD and screening paragraphs reorganized; some AMR language retained (still to remove)
- Methods: minor phrasing edits to MTCT and connector sections
- Table 1: column alignment updated; "calibrated" note removed from ANC row

All future edits should be applied to `syph_dx_zim_v2.md`.

---

## Order of operations

1. **Text-only pass (cleared, no recalibration needed):** 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 1.11, 1.14
2. **Needs debugging first:** 1.7 (investigate Fig 3 vs Fig 4 discrepancy in plot scripts)
3. **Scope decision needed:** 1.12b (budget impact analysis — decide whether to add and what to include)
4. **Recalibration batch:** 2.1 + 2.2 (early/late latent cutoff; chancre visibility + women's CS rates); then update 1.13 methods text
5. **Figure updates:** After 1.1 (remove CS from Figs 2–4) and after recalibration (update all figures)
6. **Final pass:** version number, supplement tables (Table S6 posteriors), STIsim version in methods text
