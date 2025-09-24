# 🧬 Predicting Mutation Effects with Fuzzy Logic

A fuzzy inference system for predicting the effects of genetic mutations, considering variables such as mutation type, location, population frequency, and inferring functional impact, clinical risk, and transmission probability.

# 🚀 Motivation

DNA mutations can have distinct impacts on proteins and phenotypes. Deterministic methods do not capture the nuances of these variations well, which is why a fuzzy model was chosen: it allows for dealing with uncertainty and degrees of relevance, more closely resembling human reasoning.
This project serves as an academic proof-of-concept and can be adapted to real-world bioinformatics and computational genetics scenarios.

# ⚙️ System Structure

**🔹 Inputs (antecedent linguistic variables)**

Mutation Type (e.g., synonymous, nonsynonymous, nonsense, frameshift)

Mutation Location (e.g., coding regions, introns, promoters)

Population Frequency (e.g., rare, intermediate, common)

**🔹 Outputs (consequent linguistic variables)**

Functional Impact (low, moderate, high)

Clinical Risk (low, moderate, high)

Transmission Probability (low, medium, high)

**🔹 Inference Engine**

Based on fuzzy rules (r1–r16) modeled after biomedical heuristics.
Defuzzification via the centroid method.

# 📊 Results

The system can capture nonlinear relationships between variables.
Extreme cases were tested (rare mutation + nonsense in coding region → high impact).
The results are consistent with theoretical expectations.

# 📈 Visualization

Membership functions were modeled as triangular, trapezoidal, and Gaussian.

Next steps: add interactive graphical visualization for inputs/outputs.

# 🛠️ Technologies

Python 3.11+ (https://www.python.org/)

scikit-fuzzy (https://scikit-fuzzy.readthedocs.io/en/latest/)

NumPy (https://numpy.org/)

# 📦 How to run?

```
git clone https://github.com/Icoh007/Predicting-Mutation-Effects-with-Fuzzy-Logic.git

cd fuzzy-mutation-effect-predictor

pip install -r requirements.txt

python main.py
```
