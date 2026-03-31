# Revision Plan: Syphilis Dx Zimbabwe

**Date:** 2026-03-31
**Branch:** `submission-mar31`
**Source:** Co-author feedback call (Lori Newman, Darcy Rao, Remco Peters), March 23, 2026
**Manuscript:** `syph_dx_zim_v2.md` (working file — v2 has partial edits from Lori/Darcy already applied)

---

## Status legend

- ✅ Done
- 🔄 In progress
- ⏳ Pending decision (see owner)
- ❌ Needs recalibration — do in a batch

---

## Group 1: Text-only changes (no code or recalibration needed)

| # | Change | Location | Owner | Status | Notes |
|---|--------|----------|-------|--------|-------|
| 1.1 | Remove newborn POC (CS) scenario | Abstract, Methods Table 2, Results (CS paragraph), Fig 4 caption | Robyn | ⏳ | Decided by Lori & Darcy. Removes one row from Table 2 and one results paragraph. |
| 1.2 | Standardize "pathways" → "screening and diagnostic strategies" | Throughout | Robyn | ⏳ | Check consistency: Table 1 title, section headers, results |
| 1.3 | Reframe Table 1 columns | Methods Table 1 | Robyn | ⏳ | Rename sensitivity columns to "Sensitivity (active infection)" and "False-positive rate (past/no active infection)" per feedback |
| 1.4 | Remove AMR mentions | Intro para 2; Limitations | Robyn | ⏳ | Intro: remove "possible contribution to antimicrobial resistance concerns"; Limitations: remove "nor did we model potential AMR to benzathine penicillin" |
| 1.5 | Strengthen benzathine penicillin supply angle | Policy implications | Robyn | ⏳ | Already present — move to lead the paragraph, make it more prominent |
| 1.6 | Add HIV context (1–2 sentences) | Introduction | Robyn | ⏳ | Brief mention of Zimbabwe's generalized HIV epidemic and syphilis-HIV coinfection as background |
| 1.7 | Clarify Fig 3 vs Fig 4 number discrepancy | Results: Treatment outcomes section | Robyn | ⏳ | Fig 3 shows 2000–2024 averages by pathway (~15K GUD/yr); Fig 4 shows total 2028+ (~140K total). Add a sentence explaining the time periods differ. |
| 1.8 | Update "confirmatory test" terminology | Throughout | Peters (WHO) | ⏳ | Awaiting Peters input. Suggested direction: approach-specific names, e.g. "treponemal followed by non-treponemal (T/NT) testing" |
| 1.9 | Update WHO department name for Peters affiliation | Author affiliations | Robyn | ✅ | Done in v2: "Department for HIV, TB, Hepatitis and Sexually Transmitted Infections" |
| 1.10 | Fix author middle initials (Alina M. Muellenmeister, Romesh G. Abeysuriya) | Author list | Robyn | ✅ | Done in v2 |
| 1.11 | Clarify definition of primary framing/key message | Abstract, Conclusions | Robyn | ⏳ | Awaiting Lori's steer: emphasize combined strategy + ANC bottleneck vs include undertreatment angle |
| 1.12 | Include undertreatment (brief section or expand) | Results/Discussion | Lori & all | ⏳ | Decision pending. If yes, needs modeling output; if brief, may be text-only |
| 1.13 | Care-seeking rate ~21% — add citation or revise | Methods, Results cascade | All | ⏳ | Flagged as low. Either cite source or revise estimate (would require recalibration — see Group 2) |
| 1.14 | GUD sensitivity 80% — keep + justify | Methods | All | ⏳ | Keep with explicit justification sentence; awaiting group confirmation |
| 1.15 | Role of syphilis transmission in paper — include or not | Results (Fig 1E), Discussion | Robyn | ⏳ | Currently included (Fig 1E, cascade analysis). Decision: keep brief or cut? |
| 1.16 | Newborn pathway assumptions — clarify | Methods (newborn pathway) | All | ⏳ | Clarify clinical vs surveillance definition of congenital syphilis detection |

---

## Group 2: Changes requiring recalibration

These require updating model parameters and re-running calibration + scenario analysis before the manuscript can be updated.

| # | Change | What to update | Status | Notes |
|---|--------|---------------|--------|-------|
| 2.1 | Update early/late latent definition to WHO 24 months | Currently "12–14 months" in `diseases.py`. Change cutoff from ~12 to 24 months. Affects natural history, MTCT stage distribution, screening pathway overtreatment rates | ⏳ | May shift late latent share of MTCT burden; worth checking if results change substantially before committing |
| 2.2 | Revise care-seeking rate (~21%) if decision is to change estimate | `interventions.py` care-seeking parameters | ⏳ | Only if Group 1.13 decision is to revise, not just add citation |
| 2.3 | Any parameter changes from GUD sensitivity or newborn pathway assumption decisions | `interventions.py` | ⏳ | Depends on outcomes of Group 1.14 and 1.16 |

---

## Group 3: Pending decisions (not yet assigned)

| Topic | Decision needed | Owner | Status |
|-------|----------------|-------|--------|
| Include undertreatment | Add short section vs expand modeling | Lori & all | ⏳ |
| Newborn testing scenario expansion | Add testing-at-birth scenario vs justify current | Lori & all | ⏳ |
| "Confirmatory test" renaming | Approach-specific terms per Peters suggestion | Peters | ⏳ |
| Care-seeking rate | Accept with citation vs revise | All | ⏳ |
| GUD sensitivity 80% | Keep + justify vs revert to 100% | All | ⏳ |

---

## Notes on v2 edits already made

The file `syph_dx_zim_v2.md` contains partial edits from Lori/Darcy, including:
- Title change: "novel syphilis diagnostics" → "novel tests for active syphilis"
- Abstract: minor phrasing edits, "POC" abbreviation introduced
- Intro: GUD and screening paragraphs reorganized; some AMR language retained
- Methods: minor phrasing edits to MTCT and connector sections
- Table 1: column alignment updated; "calibrated" note removed from ANC row

All future edits should be applied to `syph_dx_zim_v2.md`.

---

## Order of operations (suggested)

1. **First pass (text-only, no decisions pending):** 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7
2. **Wait for Peters input:** 1.8 (confirmatory test terminology)
3. **Wait for Lori + group decisions:** 1.11, 1.12, 1.13, 1.15, 1.16
4. **Recalibration batch:** 2.1 (early/late latent), and any others from Group 2 once decisions finalized
5. **Final pass:** update version number, figures, supplement tables after recalibration
