import csv
from collections import OrderedDict

initial_female_population = open("initial_female_population_output.csv", "r")
female_reader = csv.reader(initial_female_population, delimiter=',', quotechar='"')

female_age_collation = {}
female_population_size = 0.0

for row in female_reader:
    female_population_size += 1
    if row[1] in female_age_collation:
        female_age_collation[row[1]] += 1
    else:
        female_age_collation[row[1]] = 1

sorted_female_ages = OrderedDict(sorted(female_age_collation.items()))

print("### Female Initial Population ###")

for age in sorted_female_ages:
    print("Age: " + str(age).zfill(2)
    + ", Number: " + str(sorted_female_ages[age]).zfill(3)
    + ", Percentage: " + str(((sorted_female_ages[age]/female_population_size) * 100 )))

print("\n### Male Initial Population ###")

initial_male_population = open("initial_male_population_output.csv", "r")
male_reader = csv.reader(initial_male_population, delimiter=',', quotechar='"')

male_age_collation = {}
male_population_size = 0.0

for row in male_reader:
    male_population_size += 1
    if row[1] in male_age_collation:
        male_age_collation[row[1]] += 1
    else:
        male_age_collation[row[1]] = 1

sorted_male_ages = OrderedDict(sorted(male_age_collation.items()))

for age in sorted_male_ages:
    print("Age: " + str(age).zfill(2)
    + ", Number: " + str(sorted_male_ages[age]).zfill(3)
    + ", Percentage: " + str(((sorted_male_ages[age]/male_population_size) * 100 )))
