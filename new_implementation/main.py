import csv
import random


class Person:
    def __init__(self, age, sex, mort_dist, repr_dist, pref_dist, genes=None):
        self.age = age
        self.sex = sex
        self.alive = True
        self.mortality = mort_dist
        self.reproduction = repr_dist
        self.pref_dist = pref_dist
        if genes is None:
            self.genes = Genetics()
        else:
            self.genes = genes

    def __repr__(self):
        return ("AGE: " + str(self.age).zfill(2)
        + ", SEX: " + str(self.sex)
        + ", ALIVE: " + str(self.alive)
        + ", GENES: " + str(self.genes))

    def get_age(self):
        return self.age

    def set_age(self, age):
        self.age = age

    def increase_age(self):
        if random.uniform(0.0,1.0) < self.mortality.get_survival(self.age):
            self.age += 5
        else:
            self.alive = False

    def calculate_age_class(self):
        return ((self.age + 5) / 5)

    def set_mortality_distribution(self, distribution):
        self.mortality = distribution

    def reproduce(self):
        return self.reproduction.get_reproduction(self.age) > random.uniform(0.0,1.0)



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
        if age < len(self.probability_table) and age >= 0:
            return self.probability_table[age]
        else:
            return 0.0

    def get_survival_by_age_class(self, ageClass):
        return get_survival((ageClass * 5) - 1)

female_mortality = MortalityDistribution("female_data.csv")
#print(female_mortality)

male_mortality = MortalityDistribution("male_data.csv")
#print(male_mortality)

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
        if age < len(self.probability_table) and age >= 0:
            return self.probability_table[age]
        else:
            return 0.0

    def get_reproduction_by_age_class(self, ageClass):
        return get_reproduction((ageClass * 5) - 1)



female_reproduction = ReproductiveDistribution("female_data.csv")
#print(female_reproduction)

male_reproduction = ReproductiveDistribution("male_data.csv")
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
            self.genes = self.random_genes(8)
        else:
            self.genes = genes

    def __repr__(self):
        return str(self.genes)

    def random_genes(self, size):
        return "".join(str(random.randint(0, 1))
                for x in xrange(size))

# Set up base female population
print("### Setup Female Population ###")

female_population = []

for i in range(0,500):
    age = random.randint(0,18)*5
    female_population.append(Person(age, "f", female_mortality, female_reproduction, female_preference))

# for female in female_population:
#     print(female)
#
# print(female_population[25].reproduction)

for i in range(1,51):
    born = 0
    died = 0
    for female in female_population:
        female.increase_age()

        if not female.alive:
            female_population.remove(female)
            died += 1

        if female.reproduce():
            female_population.append(Person(0,"f", female.mortality, female.reproduction, female.pref_dist, genes=female.genes))
            born += 1

    print(str(i).zfill(2) + "  Population " + str(len(female_population)).zfill(4) + ", Born: " + str(born).zfill(3) + ", Died: " + str(died).zfill(3) + ", Change: " + str(born-died).zfill(3))

female_population.sort(key=lambda x: x.age, reverse=False)

age_groups = {}
for female in female_population:
    if female.age in age_groups:
        age_groups[female.age] += 1
    else:
        age_groups[female.age] = 1
    print(female)

for age_group in age_groups:
    print(str(age_group) + ": " +  str(age_groups[age_group]))

print("### Female Burn In Period Ended ###")

# Set up base male population
print("### Setup Male Population ###")
male_population = []

for i in range(0,500):
    age = random.randint(0,18)*5
    male_population.append(Person(age, "m", male_mortality, male_reproduction, male_map_preference))

# for male in male_population:
#     print(male)
#
# print(male_population[25].reproduction)

for i in range(1,51):
    born = 0
    died = 0
    for male in male_population:
        male.increase_age()

        if not male.alive:
            male_population.remove(male)
            died += 1

        if male.reproduce():
            male_population.append(Person(0,"m", male_mortality, male_reproduction, male_map_preference, genes=male.genes))
            born += 1

    print(str(i).zfill(2) + "  Population " + str(len(male_population)).zfill(4) + ", Born: " + str(born).zfill(3) + ", Died: " + str(died).zfill(3) + ", Change: " + str(born-died).zfill(3))

male_population.sort(key=lambda x: x.age, reverse=False)

age_groups = {}
for male in male_population:
    if male.age in age_groups:
        age_groups[male.age] += 1
    else:
        age_groups[male.age] = 1
    print(male)

for age_group in age_groups:
    print(str(age_group) + ": " +  str(age_groups[age_group]))

print("### Male Burn In Period Ended ###")

print("Female Population Size: " + str(len(female_population)))
print("Male Population Size: " + str(len(male_population)))

### Begin MAP simulation ###
print("### Begin MAP Simulation ###")
map_output = open("map_simulation_output.csv", "wb")
writer = csv.writer(map_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

map_single_females = female_population[:]
map_single_males = male_population[:]
map_couples = []

for gen in range(1,51):
    pairings = 0.0
    born = 0.0
    died = 0.0

    map_single_females.sort(key=lambda x: x.age)
    map_single_males.sort(key=lambda x: x.age)

    # Pairing Phase
    for male in map_single_males:
        for female in map_single_females:
            if female.pref_dist.get_preference(female.calculate_age_class(), male.calculate_age_class()) == 1 and male.pref_dist.get_preference(male.calculate_age_class(),female.calculate_age_class()) == 1:
                map_couples.append([female, male])
                map_single_females.remove(female)
                map_single_males.remove(male)

                pairings += 1

                break

    # Reproduction Phase
    for couple in map_couples:
        if couple[0].reproduce() + couple[1].reproduce():
            born += 1

            if 0.5 > random.uniform(0.0, 1.0):
                map_single_females.append(Person(0, "f", female_mortality, female_reproduction, female_preference, genes=couple[0].genes))
            else:
                map_single_males.append(Person(0, "m", male_mortality, male_reproduction, male_map_preference, genes=couple[1].genes))

    # Aging Phase

    for female in map_single_females:
        female.increase_age()
        if not female.alive:
            died += 1.0
            map_single_females.remove(female)

    for male in map_single_males:
        male.increase_age()
        if not male.alive:
            died += 1.0
            map_single_males.remove(male)

    for couple in map_couples:

        couple[0].increase_age()
        if not couple[0].alive:
            died += 1.0

        couple[1].increase_age()
        if not couple[1].alive:
            died += 1.0

        if not couple[0].alive and couple[1].alive:
            map_single_males.append(couple[1])
            map_couples.remove(couple)
        elif couple[0].alive and not couple[1].alive:
            map_single_females.append(couple[0])
            map_couples.remove(couple)
        elif not couple[0].alive and not couple[1].alive:
            map_couples.remove(couple)

    population = len(map_single_males) + len(map_single_females) + (len(map_couples)*2)
    print(str(gen).zfill(2) + " - Population: " + str(population).zfill(4) + ", Born: " + str(born).zfill(4) + ", Died: " + str(died).zfill(4) + ", Change: " + str(born-died).zfill(4))
