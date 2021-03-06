import csv
import random
import copy


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

    def row_builder(self):
        return [self.age, self.sex, self.alive, self.genes, self.pref_dist.name]

    def increase_age(self):
        if self.age == 85:
            self.alive = False
        elif random.uniform(0.0,1.0) < self.mortality.get_survival(self.age):
            self.age += 5
        else:
            self.alive = False

    def calculate_age_class(self):
        return ((self.age + 5) / 5)

    def reproduce(self):
        return self.reproduction.get_reproduction(self.age) > random.uniform(0.0,1.0)

    def calc_preference(self, other):
        return self.pref_dist.get_preference(self.calculate_age_class(), other.calculate_age_class())



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
            self.genes = self.random_genes(32)
        else:
            self.genes = genes

    def __repr__(self):
        return str(self.genes)

    def random_genes(self, size):
        return "".join(str(random.randint(0, 1))
                for x in xrange(size))

    def crossover(self, other):
        # Crosses over this set of genes with another set of genes with a 50% probability at each locus.
        crossed = ""

        for i in range(0,len(self.genes)):
            if 0.5 > random.uniform(0.0,1.0):
                crossed += self.genes[i]
            else:
                crossed += other.genes[i]

        return crossed

    def mutate(self, rate):
        new_genes = ""
        for gene in self.genes:
            if rate > random.uniform(0.0, 1.0):
                if gene == "0":
                    new_genes += "1"
                else:
                    new_genes += "0"
            else:
                new_genes += gene

        self.genes = new_genes

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

for i in range(1,101):
    born = 0
    died = 0
    for female in female_population:
        female.increase_age()

        if not female.alive:
            female_population.remove(female)
            died += 1

    while len(female_population) < 500:
        female = random.choice(female_population)

        if female.reproduce():
            female_population.append(Person(0,"f", female.mortality, female.reproduction, female.pref_dist, genes=female.genes))
            born += 1

    #print(str(i).zfill(2) + "  Population " + str(len(female_population)).zfill(4) + ", Born: " + str(born).zfill(3) + ", Died: " + str(died).zfill(3) + ", Change: " + str(born-died).zfill(3))

female_population.sort(key=lambda x: x.age, reverse=False) #

age_indices = []
age_groups = {}

for i in range(0,19):
    age_indices.append(i*5)

for age in age_indices:
    age_groups[age] = 0

for female in female_population:
    if female.age in age_groups:
        age_groups[female.age] += 1
    else:
        age_groups[female.age] = 1
    #print(female)

for age_group in age_indices:
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

for i in range(1,101):
    born = 0
    died = 0
    for male in male_population:
        male.increase_age()

        if not male.alive:
            male_population.remove(male)
            died += 1
    while len(male_population) < 500:
        male = random.choice(male_population)

        if male.reproduce():
            male_population.append(Person(0,"m", male_mortality, male_reproduction, male_map_preference, genes=male.genes))
            born += 1

    #print(str(i).zfill(2) + "  Population " + str(len(male_population)).zfill(4) + ", Born: " + str(born).zfill(3) + ", Died: " + str(died).zfill(3) + ", Change: " + str(born-died).zfill(3))

male_population.sort(key=lambda x: x.age, reverse=False)

age_indices = []
age_groups = {}

for i in range(0,19):
    age_indices.append(i*5)

for age in age_indices:
    age_groups[age] = 0

for male in male_population:
    if male.age in age_groups:
        age_groups[male.age] += 1
    else:
        age_groups[male.age] = 1
    #print(female)

for age_group in age_indices:
    print(str(age_group) + ": " +  str(age_groups[age_group]))

print("### Male Burn In Period Ended ###")

print("Female Population Size: " + str(len(female_population)))
print("Male Population Size: " + str(len(male_population)))

population = []

population_output = open("initial_population_output.csv", "wb")
writer = csv.writer(population_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

for female in female_population:
    population.append(female)

for male in male_population:
    population.append(male)

population.sort(key=lambda x: x.age, reverse=False)

for person in population:
    writer.writerow(person.row_builder())

age_indices = []
age_groups = {}

for i in range(0,19):
    age_indices.append(i*5)


for age in age_indices:
    age_groups[age] = 0

for person in population:
    if person.age in age_groups:
        age_groups[person.age] += 1
    else:
        age_groups[person.age] = 1

for age_group in age_indices:
    print(str(age_group) + ": " +  str(age_groups[age_group]))

### Simulation Parameters ###
mutation_rate = 0.002
simulation_generations = 1000


### Begin MAP simulation ###
print("### Begin MAP Simulation ###")
map_output = open("map_simulation_output.csv", "wb")
writer = csv.writer(map_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

map_females = copy.deepcopy(female_population)
map_males = copy.deepcopy(male_population)

output_row = [0, len(map_females), len(map_males), 0, 0, 0, 0, 0]

map_age_groups = {}

for age_index in age_indices:
    map_age_groups[age_index] = 0

for male in map_males:
    if male.age in map_age_groups:
        map_age_groups[male.age] += 1
    else:
        map_age_groups[male.age] = 1

for female in map_females:
    if female.age in map_age_groups:
        map_age_groups[female.age] += 1
    else:
        map_age_groups[female.age] = 1

for age_group in age_indices:
    output_row.append(map_age_groups[age_group])

print(output_row)
writer.writerow(output_row)

for gen in range(1,simulation_generations+1):
    born = 0.0
    died = 0.0
    female_age_reproduction = 0.0
    male_age_reproduction = 0.0
    couple_age_gap = 0.0

    oldest_female = 0
    oldest_male = 0

    map_females.sort(key=lambda x: x.age)
    map_males.sort(key=lambda x: x.age)

    # Aging Phase
    for female in map_females:
        female.increase_age()
        if not female.alive:
            died += 1.0
            map_females.remove(female)

    for male in map_males:
        male.increase_age()
        if not male.alive:
            died += 1.0
            map_males.remove(male)


    # Reproduction Phase
    while (len(map_males) + len(map_females)) < 1000:
        female = random.choice(map_females)
        male = random.choice(map_males)

        if (female.calc_preference(male) == 1 and male.calc_preference(female) == 1):
            if female.reproduce() and male.reproduce():
                # For output
                born += 1
                female_age_reproduction += female.age
                male_age_reproduction += male.age
                couple_age_gap += abs(female.age - male.age)

                if female.age > oldest_female:
                    oldest_female = female.age

                if male.age > oldest_male:
                    oldest_male = male.age

                if 0.5 > random.uniform(0.0, 1.0):
                    new_female = Person(0, "f", female_mortality, female_reproduction, female_preference)
                    new_genes = Genetics(genes=female.genes.crossover(male.genes))
                    new_genes.mutate(mutation_rate)
                    new_female.genes = new_genes
                    map_females.append(new_female)
                else:
                    new_male = Person(0, "m", male_mortality, male_reproduction, male_map_preference)
                    new_genes = Genetics(genes=male.genes.crossover(female.genes))
                    new_genes.mutate(mutation_rate)
                    new_male.genes = new_genes
                    map_males.append(new_male)

    population = len(map_males) + len(map_females)

    map_age_groups = {}

    for age_index in age_indices:
        map_age_groups[age_index] = 0

    for male in map_males:
        if male.age in map_age_groups:
            map_age_groups[male.age] += 1
        else:
            map_age_groups[male.age] = 1

    for female in map_females:
        if female.age in map_age_groups:
            map_age_groups[female.age] += 1
        else:
            map_age_groups[female.age] = 1

    female_reproductive_ave = female_age_reproduction / born
    male_reproductive_ave = male_age_reproduction / born
    ave_age_gap = couple_age_gap / born

    output_string = "GEN: " + str(gen).zfill(2)
    output_string += ", Females: " + str(len(map_females))
    output_string += ", Males: " + str(len(map_males))
    output_string += ", Fem Rep Age: " + str(female_reproductive_ave)
    output_string += ", Mal Rep Age: " + str(male_reproductive_ave)
    output_string += ", Couple Age Gap: " + str(ave_age_gap)
    output_string += ", Oldest Female Rep: " + str(oldest_female)
    output_string += ", Oldest Male Rep: " + str(oldest_male) + ", "

    output_row = [gen, len(map_females), len(map_males), female_reproductive_ave, male_reproductive_ave, ave_age_gap, oldest_female, oldest_male]

    for age_group in age_indices:
        output_string += str(age_group).zfill(2) + ": " +  str(map_age_groups[age_group]).zfill(3) + ", "
        output_row.append(map_age_groups[age_group])
    #print(output_string)

    writer.writerow(output_row)

map_population = []

map_population_output = open("map_population_output.csv", "wb")
writer = csv.writer(map_population_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

for female in map_females:
    map_population.append(female)

for male in map_males:
    map_population.append(male)

map_population.sort(key=lambda x: x.age, reverse=False)

for person in map_population:
    writer.writerow(person.row_builder())

### Begin MYP simulation ###
print("### Begin MYP Simulation ###")
myp_output = open("myp_simulation_output.csv", "wb")
writer = csv.writer(myp_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

myp_females = copy.deepcopy(female_population)
myp_males = copy.deepcopy(male_population)

for male in myp_males:
    male.pref_dist = male_myp_preference

output_row = [0, len(myp_females), len(myp_males), 0, 0, 0, 0, 0]

myp_age_groups = {}

for age_index in age_indices:
    myp_age_groups[age_index] = 0

for male in myp_males:
    if male.age in myp_age_groups:
        myp_age_groups[male.age] += 1
    else:
        myp_age_groups[male.age] = 1

for female in myp_females:
    if female.age in myp_age_groups:
        myp_age_groups[female.age] += 1
    else:
        myp_age_groups[female.age] = 1

for age_group in age_indices:
    output_row.append(myp_age_groups[age_group])

print(output_row)
writer.writerow(output_row)

for gen in range(1,simulation_generations+1):
    born = 0.0
    died = 0.0
    female_age_reproduction = 0.0
    male_age_reproduction = 0.0
    couple_age_gap = 0.0

    oldest_female = 0
    oldest_male = 0

    myp_females.sort(key=lambda x: x.age)
    myp_males.sort(key=lambda x: x.age)

    # Aging Phase
    for female in myp_females:
        female.increase_age()
        if not female.alive:
            died += 1.0
            myp_females.remove(female)

    for male in myp_males:
        male.increase_age()
        if not male.alive:
            died += 1.0
            myp_males.remove(male)


    # Reproduction Phase
    while (len(myp_males) + len(myp_females)) < 1000:
        female = random.choice(myp_females)
        male = random.choice(myp_males)

        if (female.calc_preference(male) == 1 and male.calc_preference(female) == 1):
            if female.reproduce() and male.reproduce():
                # For output
                born += 1
                female_age_reproduction += female.age
                male_age_reproduction += male.age
                couple_age_gap += abs(female.age - male.age)

                if female.age > oldest_female:
                    oldest_female = female.age

                if male.age > oldest_male:
                    oldest_male = male.age

                if 0.5 > random.uniform(0.0, 1.0):
                    new_female = Person(0, "f", female_mortality, female_reproduction, female_preference)
                    new_genes = Genetics(genes=female.genes.crossover(male.genes))
                    new_genes.mutate(mutation_rate)
                    new_female.genes = new_genes
                    myp_females.append(new_female)
                else:
                    new_male = Person(0, "m", male_mortality, male_reproduction, male_myp_preference)
                    new_genes = Genetics(genes=male.genes.crossover(female.genes))
                    new_genes.mutate(mutation_rate)
                    new_male.genes = new_genes
                    myp_males.append(new_male)

    population = len(myp_males) + len(myp_females)

    myp_age_groups = {}

    for age_index in age_indices:
        myp_age_groups[age_index] = 0

    for male in myp_males:
        if male.age in myp_age_groups:
            myp_age_groups[male.age] += 1
        else:
            myp_age_groups[male.age] = 1

    for female in myp_females:
        if female.age in myp_age_groups:
            myp_age_groups[female.age] += 1
        else:
            myp_age_groups[female.age] = 1

    female_reproductive_ave = female_age_reproduction / born
    male_reproductive_ave = male_age_reproduction / born
    ave_age_gap = couple_age_gap / born

    output_string = "GEN: " + str(gen).zfill(2)
    output_string += ", Females: " + str(len(myp_females))
    output_string += ", Males: " + str(len(myp_males))
    output_string += ", Fem Rep Age: " + str(female_reproductive_ave)
    output_string += ", Mal Rep Age: " + str(male_reproductive_ave)
    output_string += ", Couple Age Gap: " + str(ave_age_gap)
    output_string += ", Oldest Female Rep: " + str(oldest_female)
    output_string += ", Oldest Male Rep: " + str(oldest_male) + ", "

    output_row = [gen, len(myp_females), len(myp_males), female_reproductive_ave, male_reproductive_ave, ave_age_gap, oldest_female, oldest_male]

    for age_group in age_indices:
        output_string += str(age_group).zfill(2) + ": " +  str(myp_age_groups[age_group]).zfill(3) + ", "
        output_row.append(myp_age_groups[age_group])

    #print(output_string)

    writer.writerow(output_row)

myp_population = []

myp_population_output = open("myp_population_output.csv", "wb")
writer = csv.writer(myp_population_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

for female in myp_females:
    myp_population.append(female)

for myle in myp_males:
    myp_population.append(male)

myp_population.sort(key=lambda x: x.age, reverse=False)

for person in myp_population:
    writer.writerow(person.row_builder())
