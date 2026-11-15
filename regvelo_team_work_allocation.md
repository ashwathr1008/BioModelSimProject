# RegVelo Project: Five-Person Work Allocation

**Project:** How does incomplete developmental sampling affect RegVelo's predictions of cell fate and regulatory drivers in human muscle development?

**Duration:** Two weeks  
**Team:** Five people  
**Planned core experiments:** Three conditions (full data, targeted intermediate-cell removal, and matched random removal), three repeat blocks, nine model fits, and two separate regulator perturbations per fit (MSC and MYOG).

Assign ownership by pipeline component. Each person develops one shared implementation that handles every experimental condition. This keeps the comparison consistent and lets downstream analysis proceed as individual models finish.

## 1. Responsibilities

| Person | Main responsibility | Work that can start immediately | Main deliverables |
|---|---|---|---|
| **1 — Data, sampling, and preprocessing** | Inspect inputs, construct sampling conditions, and prepare model inputs. | Verify spliced/unspliced input provenance; inspect population, stage, and sample counts; implement targeted and random removal; check gene and regulatory-prior alignment. | Prepared datasets, retained/removed cell-ID manifests, preprocessing configuration, and sampling/QC figures. |
| **2 — Environment, baseline, and model fitting** | Establish the environment, reproduce the tutorial, maintain the shared training script, and coordinate model runs. | Run the existing tutorial while Person 1 develops experimental preprocessing; record package versions, training settings, runtime, and memory use. | Reproduced baseline, reusable fitting script, saved models, velocities, latent times, and run logs. |
| **3 — Cell-fate and trajectory analysis** | Implement shared fate inference and trajectory measurements. | Inspect annotated transitions and endpoints; develop the terminal-anchor selection procedure; write fate-inference and trajectory-scoring functions. | Fixed terminal-anchor manifest, shared fate-inference function, baseline fate probabilities, and trajectory results. |
| **4 — MSC/MYOG perturbation analysis** | Implement and run the two regulator perturbations for every fitted model. | Build the perturbation wrapper using the tutorial; verify the meaning of its outputs; document biological expectations and coordinate interfaces with Person 3. | Perturbed fate probabilities and per-cell changes relative to the matching unperturbed model. |
| **5 — Comparisons, statistics, and figures** | Align comparable cells, calculate condition contrasts, aggregate repeats, and produce consistent figures. | Define result-table schemas; implement alignment and comparison functions using illustrative inputs; prepare plotting functions and checks. | Combined results tables, targeted-versus-random comparisons, repeat-level effect-size plots, and figure-generation scripts. |

Person 5 has a substantial analysis role from day one. Writing and presentation preparation are shared by all five people.

## 2. Component Boundaries and Handoffs

### Person 1 → Person 2: model inputs

Provide each dataset with:

- Condition and repeat identifiers.
- Exact retained and removed cell IDs.
- Sampling seed, population counts, and removal count.
- Gene order and the aligned regulatory prior.
- Preprocessing settings and checks that the inputs match the agreed project plan.

Apply the same preprocessing procedure independently to each condition. Rebuild the required neighborhood-dependent quantities after cell removal.

### Person 2 → Persons 3 and 4: fitted models

Release each completed fit immediately with:

- Model checkpoint and the corresponding prepared dataset.
- Velocity and latent-time outputs.
- Training configuration, model seed, package versions, and completion status.
- Runtime and any warnings or failed checks.

Persons 3 and 4 should not wait for all nine fits before starting real analyses.

### Person 3 ↔ Person 4: shared fate inference

Person 3 owns the function used to infer fates from both unperturbed and perturbed outputs. Person 4 calls that function rather than maintaining a separate version.

Choose terminal anchors once using the agreed baseline procedure, save their exact cell IDs, and protect the required endpoint populations during sampling. Keep endpoint names, ordering, and fate-inference settings consistent across runs.

### Persons 3 and 4 → Person 5: analysis outputs

Provide cell-level tables indexed by condition, repeat, and cell ID, with explicit endpoint and regulator labels where applicable.

Person 5:

- Restricts the primary comparison to common retained nonterminal cells within each repeat block.
- Calculates each regulator's effect relative to its own fitted model's unperturbed prediction.
- Compares targeted and random removal against the corresponding full-data reference.
- Summarizes variation across repeat blocks without treating computational repeats or individual cells as independent biological replicates.
- Reports missing outputs and failed runs explicitly.

Persons 3–5 can develop their functions before the main fits exist. Illustrative inputs are for checking code; all reported scientific results must come from actual experimental runs.

## 3. Two-Week Schedule

| Period | Parallel work | Completion target |
|---|---|---|
| **Days 1–2** | Person 1 inspects and prepares inputs. Person 2 establishes the environment and runs the tutorial pilot. Persons 3–5 develop downstream functions, interfaces, and checks. | Working tutorial pilot, agreed experiment configuration, and compatible component interfaces. |
| **Days 3–4** | Connect the components and complete the first full/targeted/random repeat block, including both perturbations and comparison plots. All owners resolve integration problems in their components. | One complete repeat block processed from input data through final comparisons. |
| **Days 4–8** | Person 2 coordinates remaining fits. Persons 3–5 process completed runs as they arrive. Person 1 checks sampling consistency and helps diagnose input-related failures. | Planned model fits and most downstream analyses complete. |
| **Days 9–10** | Resolve failed or missing outputs, finalize comparison tables and figures, and interpret results together. | Agreed results and limitations; core analysis frozen. |
| **Days 11–14** | Each person writes their methods/results contribution and prepares corresponding slides. Assemble, review, and rehearse together. | Final report, presentation, and reproducible project package. |

The first complete repeat block is the integration milestone. Validate it before launching the remaining main experiments.

## 4. Coordination Rules

1. **Agree on inputs and outputs on day one.** Use consistent condition names, repeat IDs, cell IDs, endpoint names, and regulator labels. Preserve model and preprocessing settings with every run.
2. **Use one shared implementation across conditions.** Do not assign separate people to build independent full-data, targeted-removal, and random-removal pipelines.
3. **Separate sampling seeds from model seeds.** Record both, and use the agreed repeat-block pairing so that comparisons can be reconstructed.
4. **Hand off runs individually.** Training, fate inference, perturbation analysis, and aggregation should overlap whenever their dependencies allow.
5. **Coordinate compute separately from ownership.** With one GPU, Person 2 maintains a common queue for training and GPU-dependent perturbation jobs. With multiple GPUs, distribute identical jobs using the same pinned environment and configuration.
6. **Assign a backup for each component.** The backup should be able to run or resume the component from its instructions and saved outputs.
7. **Keep a shared run manifest.** Track planned, running, complete, failed, and analyzed jobs, together with output locations and brief failure notes.
8. **Keep the two-week scope fixed.** Complete the core experiment before considering another tissue, extra regulators, or additional removal levels.

## 5. AI Assistance

Each owner can use AI to draft code, explain unfamiliar APIs, debug errors, and prepare documentation within their component.

The owner remains responsible for verifying API behavior, input/output alignment, and scientific interpretation. Review generated code before expensive runs. Use the first completed repeat block to check that the components agree.

## 6. Writing and Presentation Ownership

| Person | Primary writing and slide contribution |
|---|---|
| **1** | Dataset, sampling design, preprocessing, and sampling/QC figures. |
| **2** | RegVelo setup, baseline reproduction, model fitting, and computational feasibility. |
| **3** | Trajectory and fate methods, results, and interpretation. |
| **4** | MSC/MYOG perturbation methods, results, and biological context. |
| **5** | Comparison methods, repeat-level variability, consolidated figures, and consistency of reported numbers. |

Everyone contributes to the introduction, discussion, limitations, and rehearsal. Choose one person to assemble the final report and one to assemble the slides, while each component owner supplies finished text and figures.

## 7. Completion Checklist

- [ ] Tutorial reproduction and environment documented.
- [ ] Sampling conditions and retained/removed cell IDs saved.
- [ ] Planned fits completed, or failures and scope limitations documented.
- [ ] Shared fate inference used consistently for unperturbed and perturbed outputs.
- [ ] MSC and MYOG analyses completed for each usable fitted model.
- [ ] Primary comparisons use the agreed common-cell evaluation sets.
- [ ] Figures and numerical results can be regenerated from saved outputs.
- [ ] Report distinguishes stability relative to the full-data model from biological accuracy.
- [ ] All five members review the final report and presentation.

