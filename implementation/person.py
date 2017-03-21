import csv
import random

# Import female life data table
female_life_data = []

with open("female_data.csv", "rb") as f:
    reader = csv.reader(f)
    for row in reader:
        female_life_data.append({"survival": float(row[3]), "fertility": float(row[2])})

# for row in female_life_data:
#     print(row)

# Import male life data table
male_life_data = []

with open("female_data.csv", "rb") as f:
    reader = csv.reader(f)
    for row in reader:
        male_life_data.append({"survival": float(row[3]), "fertility": float(row[2])})

# Import female preference
female_preference = []

with open("female_pref.csv", "rb") as f:
    reader = csv.reader(f)
    for row in reader:
        new_row = []
        for cell in row:
            new_row.append(int(cell))
        female_preference.append(new_row)

# Import male preferences
male_map = []

with open("male_map.csv", "rb") as f:
    reader = csv.reader(f)
    for row in reader:
        new_row = []
        for cell in row:
            new_row.append(int(cell))
        male_map.append(new_row)

male_myp = []

with open("male_myp.csv", "rb") as f:
    reader = csv.reader(f)
    for row in reader:
        new_row = []
        for cell in row:
            new_row.append(int(cell))
        male_myp.append(new_row)

class Person:
    def __init__(self, sex, age, genetics):
        # initialise class
        self.sex = sex
        self.age = age
        self.genetics = genetics
        self.alive = True

        self.fitness = 0.0
        self.calculate_fitness()

        self.mating_preferences = []
        self.calculate_mating_preferences()


    def __repr__(self):
        # provide printable representation of the class
        return "Sex: " + self.sex + ", Age: " + str(self.age) + ", Genetics: " + self.genetics + ", Fitness: " + str(self.fitness) + ", Alive: " + str(self.alive)

    def get_list(self):
        # returns a list for csv output
        return [self.sex, self.age, str(self.alive)]

    def calculate_fitness(self):
        # calculates fitness of individual
        pass

    def age_person(self):
        # looks up age in life table and either kills individual or ages them up
        life_chance = 0.0
        if self.sex == "f":
            life_chance = female_life_data[self.age]["survival"]
        elif self.sex == "m":
            life_chance = male_life_data[self.age]["survival"]
        else:
            life_chance = female_lif_data[self.age]["survival"]
        if random.uniform(0.0,1.0) < life_chance and self.alive == True:
            self.age += 5
        else:
            self.alive = False

    def crossover(self, other):
        # crosses over genetic information between this individual and another for producing offspring.
        pass

    def mutate(self):
        # mutates genetic string_length
        pass

    def reproduce(self):
        if self.sex == "f":
            return female_life_data[self.age]["fertility"]
        else:
            return female_life_data[self.age]["fertility"]
#
    def calculate_mating_preferences(self):
        if self.sex == "f":
            self.mating_preferences = female_preference
        elif self.sex == "m":
            self.mating_preferences = male_map
        else:
            self.mating_preferences = male_myp

    def set_mating_preference(self, pref):
        if self.sex == "f" or pref.lower() == "female":
            self.mating_preferences = female_preference
        elif pref.lower() == "map":
            self.mating_preferences = male_map
        elif pref.lower() == "myp":
            self.mating_preferences = male_myp
        else:
            self.mating_preferences = female_preference

    def calculate_age_class(self):
        return ((self.age + 5) / 5)

    def get_preference(self, other):
        return self.mating_preferences[self.calculate_age_class()][other.calculate_age_class()]



#
# joanne = Person("f", 22, "00101110")
#
# print(joanne)
#
# joanne.genetics = "1100000"
#
# joanne.calculate_fitness()
#
# print (joanne)
