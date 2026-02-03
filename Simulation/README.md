# Simulation Data Generation

This document describes how we generated simulated transcriptomic data to validate the model.

- Simulator: beers (v1.0)
- Transcript annotations: human genome GTF (to build the transcriptome space)
- Full configuration: see Simulation/config.yaml

We evaluated the impact of common biases by adjusting simulator parameters. The following settings were used in bias experiments:
**1) Mapping errors**
- Introduce base substitutions/indels to increase alignment difficulty.

```yaml
substitution_rate: 0.005
deletion_rate: 0.0001
insertion_rate: 0.0001
```

**2) 5'–3' positional bias**
- Simulate coverage bias along transcript length (first-strand); reduced priming density.

```yaml
primes_per_kb: 30  # reduced from the default 50
```

**3) PCR / GC bias**
- Simulate amplification bias related to GC content.

```yaml
gc_bias_constant: 0.8       # reduced base value
gc_bias_linear: 0.02        # linear term
gc_bias_quadratic: -0.01    # negative => downward-opening quadratic
```

Notes
- The defaults in Simulation/config.yaml currently show un-biased values (e.g., primes_per_kb: 50 and error rates set to 0.0). Modify those fields to the values above to reproduce the bias experiments.
- Key locations in the config:
	- Positional bias: first_strand_synthesis_step.primes_per_kb
	- Mapping errors: pcr_amplification_step.[substitution_rate|deletion_rate|insertion_rate]
	- GC bias: pcr_amplification_step.[gc_bias_constant|gc_bias_linear|gc_bias_quadratic]

