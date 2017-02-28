import csv
import random


class Person:
    def __init__(self, age, sex, mort_dist, fert_dist, pref_dist, genes=None):
        self.age = age
        self.sex = sex
        self.alive = true
        self.mortality = mort_dist
        self.fertility = fert_dist
        self.pref_dist = pref_dist
        if genes is None:
            self.genes = Genetics()
        else:
            self.genes = genes

    def __repr__(self):
        return "AGE: " + str(self.age).zfill(2) + ", SEX: " + str(self.sex) + ", ALIVE: " + str(self.alive)

    def get_age(self):
        return self.age

    def set_age(self, age):
        self.age = age

    def increase_age(self):
        if random.uniform(0.0,1.0) < self.mortality.get_survival(self.age):
            self.age += 5
        else:
            self.alive = false

    def calculate_age_class(self):
        return ((self.age + 5) / 5)

    def set_mortality_distribution(self, distribution):
        self.mortality = distribution



# Set up MortalityDistribution

class MortalityDistribution:
    def __init__(self, input_file):
        self.probability_table = []
        with open(input_file, "rb") as f:
            reader = csv.reader(f)
            for row in reader:
                self.probability_table.append(float(row[3]))

    def __repr__(self):
        str_representation = ""
        for row in self.probability_table:
            str_representation += str(row).zfill(0) + "\n"
        return str_representation

    def get_survival(self, age):
        if age < len(probability_table) and age >= 0:
            return probability_table[age]
        else:
            return 0.0

    def get_survival_by_age_class(self, ageClass):
        return get_survival((ageClass * 5) - 1)

female_fertility = MortalityDistribution("female_data.csv")
#print(female_fertility)

male_fertility = MortalityDistribution("male_data.csv")
#print(male_fertility)

print("Initial mortality distributions set up")

# Set up ReproductiveDistribution

class ReproductiveDistribution:
    def __init__(self, input_file):
        self.probability_table = []
        with open(input_file, "rb") as f:
            reader = csv.reader(f)
            for row in reader:
                self.probability_table.append(float(row[2]))

    def __repr__(self):
        str_representation = ""
        for row in self.probability_table:
            str_representation += str(row).zfill(0) + "\n"
        return str_representation

    def get_reproduction(self, age):
        if age < len(probability_table) and age >= 0:
            return probability_table[age]
        else:
            return 0.0

    def get_reproduction_by_age_class(self, ageClass):
        return get_reproduction((ageClass * 5) - 1)



female_reproduction = MortalityDistribution("female_data.csv")
#print(female_reproduction)

male_reproduction = MortalityDistribution("male_data.csv")
#print(male_reproduction)
print("Initial reproductive distributions set up")

# Set up PreferenceDistribution
class PreferenceDistribution:
    def __init__(self, name, input_file):
        self.name = name
        self.distribution = []
        with open(input_file, "rb") as f:
            reader = csv.reader(f)
            for row in reader:
                new_row = []
                for cell in row:
                    new_row.append(int(cell))
                self.distribution.append(new_row)

    def __repr__(self):
        str_representation = ""

        for row in self.distribution:
            row_str = ""
            for cell in row:
                row_str += str(cell) + ", "
            row_str += "\n"

            str_representation += row_str

        return str_representation

    def get_preference(self, age_class_one, age_class_two):
        return self.distribution[age_class_one][age_class_two]

female_preference = PreferenceDistribution("female", "female_pref.csv")
#print(female_preference)

male_map_preference = PreferenceDistribution("male_map", "male_map.csv")
#print(male_map_preference)

male_myp_preference = PreferenceDistribution("male_myp", "male_myp.csv")
#print(male_myp_preference)

print("Initial preference distributions set up")

class Genetics:
    def __init__(self, genes=None):
        if genes is None:
            self.genes = random_genes(8)
        else:
            self.genes = genes

    def random_genes(self, size):
        return "".join(str(random.randint(0, 1))
                for x in xrange(size))
