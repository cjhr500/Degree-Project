import csv
import random


class Person:
    def __init__(self, age, sex):
        self.age = age
        self.sex = sex
        self.alive = true

    def __repr__(self):
        return "AGE: " + str(self.age).zfill(2) + ", SEX: " + str(self.sex) + ", ALIVE: " + str(self.alive)

    def get_age(self):
        return self.age

    def set_age(self, age):
        self.age = age

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

print("Mortality distribution set up")

female_fertility = MortalityDistribution("female_data.csv")
print(female_fertility)

male_fertility = MortalityDistribution("male_data.csv")
print(male_fertility)

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

print("Reproductive distribution set up")

female_reproduction = MortalityDistribution("female_data.csv")
print(female_reproduction)

male_reproduction = MortalityDistribution("male_data.csv")
print(male_reproduction)
