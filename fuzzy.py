# @title Predicting the Effects of Genetic Mutations
import numpy as np
import skfuzzy as fuzzy
from skfuzzy import control as ctrl
import pandas as pd
import matplotlib.pyplot as plt
import os

# @title Antecedents
tip_mut = ctrl.Antecedent(np.arange(0, 101, 0.5), 'Gene Mutations Types')
loc_mut = ctrl.Antecedent(np.arange(0, 101, 0.5), 'Mutation Location')
freq_pop = ctrl.Antecedent(np.arange(0, 101, 0.5), 'Population Frequency')

# @title Consequents
impact = ctrl.Consequent(np.arange(0, 101, 0.5), 'Expected Functional Impact')
risk = ctrl.Consequent(np.arange(0, 101, 0.5), 'Associated Clinical Risk')
probability = ctrl.Consequent(
    np.arange(0, 101, 0.5), 'Probability of Relevant Transmission')

# @title Universe of Gene Mutations Types
tip_mut['SNP'] = fuzzy.trapmf(tip_mut.universe, [0, 0, 10, 30])
tip_mut['Insertion'] = fuzzy.trapmf(tip_mut.universe, [20, 35, 45, 60])
tip_mut['Deletion'] = fuzzy.trapmf(tip_mut.universe, [55, 75, 100, 100])

# @title Universe of Mutation Location
loc_mut['Non-coding'] = fuzzy.trapmf(loc_mut.universe, [0, 0, 10, 25])
loc_mut['Promoter'] = fuzzy.trapmf(loc_mut.universe, [15, 30, 40, 55])
loc_mut['Regulatory'] = fuzzy.trapmf(loc_mut.universe, [45, 55, 65, 75])
loc_mut['Coding_critical'] = fuzzy.trapmf(loc_mut.universe, [70, 85, 100, 100])

# @title Universe of Population Frequency
freq_pop['Very common'] = fuzzy.gaussmf(freq_pop.universe, 15, 12)
freq_pop['Common'] = fuzzy.gaussmf(freq_pop.universe, 35, 12)
freq_pop['Moderate'] = fuzzy.gaussmf(freq_pop.universe, 50, 12)
freq_pop['Rare'] = fuzzy.gaussmf(freq_pop.universe, 75, 8)
freq_pop['Very rare'] = fuzzy.gaussmf(freq_pop.universe, 95, 8)

# @title Universe of Expected Functional Impact
impact['Very low'] = fuzzy.gaussmf(impact.universe, 10, 8)
impact['Low'] = fuzzy.gaussmf(impact.universe, 30, 8)
impact['Moderate'] = fuzzy.gaussmf(impact.universe, 50, 10)
impact['High'] = fuzzy.gaussmf(impact.universe, 65, 12)
impact['Very high'] = fuzzy.gaussmf(impact.universe, 80, 14)

# @title Universe of Associated Clinical Risk
risk['Neutral'] = fuzzy.trapmf(risk.universe, [0, 0, 10, 30])
risk['Predisposition'] = fuzzy.trapmf(risk.universe, [20, 45, 60, 75])
risk['Pathogenic'] = fuzzy.trapmf(risk.universe, [60, 75, 100, 100])

# @title Universe of Relevant Transmission Probability
probability['Very low'] = fuzzy.gaussmf(probability.universe, 10, 8)
probability['Low'] = fuzzy.gaussmf(probability.universe, 30, 10)
probability['Medium'] = fuzzy.gaussmf(probability.universe, 50, 12)
probability['High'] = fuzzy.gaussmf(probability.universe, 70, 10)
probability['Very high'] = fuzzy.gaussmf(probability.universe, 90, 8)

# @title Plotting
tip_mut.view()
loc_mut.view()
freq_pop.view()

impact.view()
risk.view()
probability.view()

# @title Defining Inferences
r1 = ctrl.Rule(tip_mut['Deletion'] & loc_mut['Coding_critical'] & freq_pop['Very rare'],
               (impact['Very high'], risk['Pathogenic'], probability['High']))

r2 = ctrl.Rule(tip_mut['Insertion'] & loc_mut['Coding_critical'] & freq_pop['Very rare'],
               (impact['Very high'], risk['Pathogenic'], probability['High']))

r3 = ctrl.Rule(tip_mut['SNP'] & loc_mut['Coding_critical'] & freq_pop['Very rare'],
               (impact['High'], risk['Predisposition'], probability['Medium']))

r4 = ctrl.Rule(tip_mut['SNP'] & loc_mut['Coding_critical'] &
               freq_pop['Very common'], (impact['Low'], risk['Neutral'], probability['Low']))

r5 = ctrl.Rule(tip_mut['Deletion'] & loc_mut['Regulatory'] & freq_pop['Rare'],
               (impact['Moderate'], risk['Predisposition'], probability['Medium']))

r6 = ctrl.Rule(tip_mut['Insertion'] & loc_mut['Promoter'] & freq_pop['Rare'],
               (impact['High'], risk['Predisposition'], probability['Medium']))

r7 = ctrl.Rule(tip_mut['SNP'] & loc_mut['Promoter'] & freq_pop['Moderate'],
               (impact['Moderate'], risk['Predisposition'], probability['Medium']))

r8 = ctrl.Rule(loc_mut['Non-coding'] & freq_pop['Very common'],
               (impact['Very low'], risk['Neutral'], probability['Low']))

r9 = ctrl.Rule(loc_mut['Non-coding'] & freq_pop['Very rare'],
               (impact['Low'], risk['Predisposition'], probability['Low']))

r10 = ctrl.Rule(tip_mut['Insertion'] & loc_mut['Regulatory'] &
                freq_pop['Common'], (impact['Low'], risk['Neutral'], probability['Low']))

r11 = ctrl.Rule(tip_mut['Deletion'] & loc_mut['Promoter'] & freq_pop['Moderate'],
                (impact['High'], risk['Predisposition'], probability['Medium']))

r12 = ctrl.Rule(tip_mut['SNP'] & loc_mut['Regulatory'] & freq_pop['Rare'],
                (impact['Moderate'], risk['Predisposition'], probability['Low']))

r13 = ctrl.Rule(tip_mut['SNP'] & loc_mut['Promoter'] & freq_pop['Very common'],
                (impact['Low'], risk['Neutral'], probability['Low']))

r14 = ctrl.Rule(tip_mut['Insertion'] & loc_mut['Coding_critical'] & freq_pop['Common'],
                (impact['Moderate'], risk['Predisposition'], probability['Medium']))

r15 = ctrl.Rule(tip_mut['Deletion'] & loc_mut['Coding_critical'] & freq_pop['Common'],
                (impact['High'], risk['Predisposition'], probability['Medium']))

r16 = ctrl.Rule((tip_mut['Insertion'] | tip_mut['Deletion']) & (loc_mut['Coding_critical'] | loc_mut['Promoter'])
                & freq_pop['Very rare'], (impact['Very high'], risk['Pathogenic'], probability['High']))

# @title Creating a "rules" array
rules = [globals()[f"r{i}"] for i in range(1, 17)]

# @title Creating the Inference Engine
system = ctrl.ControlSystem(rules)


def run_case(tip_val, loc_val, fre_val):
    sim = ctrl.ControlSystemSimulation(system)

    sim.input['Gene Mutations Types'] = tip_val
    sim.input['Mutation Location'] = loc_val
    sim.input['Population Frequency'] = fre_val

    sim.compute()

    impact.view(sim=sim)
    risk.view(sim=sim)
    probability.view(sim=sim)

    out_impact = sim.output['Expected Functional Impact']
    out_risk = sim.output['Associated Clinical Risk']
    out_probability = sim.output['Probability of Relevant Transmission']

    print(
        f"INPUTS:\n\tGene Mutations Types: {tip_val:.2f}\n\n\tMutation Location: {loc_val:.2f}\n\n\tPopulation Frequency: {fre_val:.2f}\n\n")
    print(
        f"OUTPUTS:\n\tExpected Functional Impact: {out_impact:.2f}\n\n\tAssociated Clinical Risk: {out_risk:.2f}\n\n\tProbability of Relevant Transmission: {out_probability:.2f}\n")


# @title Batch testing
batch_testing = [
    (80, 85, 90)
    (89, 6, 25),
    (18, 70, 53),
    (1, 92, 60),
    (41, 62, 87),
    (42, 41, 64),
    (64, 43, 5),
    (16, 25, 53),
    (74, 65, 81),
    (74, 24, 59),
    (79, 76, 80),
    (24, 25, 88),
    (40, 25, 22),
    (13, 57, 44),
    (12, 8, 72),
    (63, 74, 82),
    (32, 84, 96),
    (7, 28, 43),
    (0, 48, 27),
    (54, 56, 24),
    (84, 70, 88)
]

for i, (mut_type, mut_loc, pop_freq) in enumerate(batch_testing, 1):
    run_case(mut_type, mut_loc, pop_freq)
    print(f"--- Caso {i} concluído ---\n")
