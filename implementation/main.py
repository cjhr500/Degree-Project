from person import Person
import random
import csv

def random_binary_string(size):
    # generates random binary string of length = size
    return "".join(str(random.randint(0, 1))
                for x in xrange(size))

### Setup Female Population ###
print("### Setup Female Population ###")

female_population = []
old_female_population_size = 0

for i in range(0,500):
    female_population.append(Person("f", 0, random_binary_string(8)))
    old_female_population_size += 1

# born = 0.0
# age_of_reproduction = 0.0
# oldest_reproduction = 0.0
#
# died = 0.0
# age_at_death = 0.0
# oldest_death = 0.0
# population_max = 0.0

total_females_born = 0
total_females_died = 0

for i in range(1,51):
    # if len(population) > population_max:
    #     population_max = len(population)
    born = 0
    died = 0
    for person in female_population:

        person.age_person()

        if not person.alive:
            # age_at_death += person.age
            #
            # if person.age > oldest_death:
            #     oldest_death = person.age
            #
            died += 1
            female_population.remove(person)

        if person.reproduce() > random.uniform(0.0, 1.0):
            # age_of_reproduction += person.age
            #
            # if person.age > oldest_reproduction:
            #     oldest_reproduction = person.age
            female_population.append(Person("f", 0, person.genetics))
            born += 1


    print("Iteration: " + str(i).zfill(3) + ", Population size: " + str(len(female_population)).zfill(4) + ", Born: " + str(born).zfill(3) + ", Died: " + str(died).zfill(3) + ", Change: " + str(len(female_population) - old_female_population_size).zfill(4))

    total_females_born += born
    total_females_died += died

    old_female_population_size = len(female_population)

    if len(female_population) == 0:
        break

print("Total born: " + str(total_females_born) + ", Total died: " + str(total_females_died))
    # if len(population) > population_max:
    #     population_max = len(population)
female_population.sort(key=lambda x: x.age, reverse=False) # sort population by age

initial_female_population = open("initial_female_population_output.csv", "wb")
writer = csv.writer(initial_female_population, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

for female in female_population:
    writer.writerow(female.get_string_list())

initial_female_population.close()
# ages = []
#
# for person in population:
#     print(person)
#     ages.append(person.age)
#
# ages.sort()
#
# print ages
#
# av_age_of_reproduction = age_of_reproduction / born
# av_age_at_death = age_at_death / died

# print("Born: " + str(born) + ", Average age of reproduction: " + str(av_age_of_reproduction) + ", Oldest reproduction: " + str(oldest_reproduction))
# print("Died: " + str(died) + ", Average age at death: " + str(av_age_at_death) + ", Oldest death: " + str(oldest_death))
# print("Max population size reached: " + str(population_max) + ", Alive at end: " + str(len(population)))

print("### Female Burn In Period Ended ###")

### Setup male population ###
print("### Setup Male Population ###")

male_population = []
old_male_population_size = 0

for i in range(0,500):
    male_population.append(Person("m", 0, random_binary_string(8)))
    old_male_population_size += 1

total_males_born = 0
total_males_died = 0

for i in range(1,51):
    born = 0
    died = 0
    for person in male_population:
        person.age_person()

        if not person.alive:
            died += 1
            male_population.remove(person)

        if person.reproduce() > random.uniform(0.0, 1.0):
            born += 1
            male_population.append(Person("m", 0, person.genetics))

    print("Iteration: " + str(i).zfill(3) + ", Population size: " + str(len(male_population)).zfill(4) + ", Born: " + str(born).zfill(3) + ", Died: " + str(died).zfill(3) + ", Change: " + str(len(male_population) - old_male_population_size).zfill(4))

    total_males_born += born
    total_males_died += died

    old_male_population_size = len(male_population)

    if len(male_population) == 0:
        break

print("Total born: " + str(total_males_born) + ", Total died: " + str(total_males_died))

male_population.sort(key=lambda x: x.age, reverse=False) # sort population by age

initial_male_population = open("initial_male_population_output.csv", "wb")
writer = csv.writer(initial_male_population, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

for male in male_population:
    writer.writerow(male.get_string_list())

initial_female_population.close()

print("### Male Burn In Period Ended ###")

### Begin MAP simulation ###
print("### Begin MAP Simulation ###")

map_output = open("map_simulation_output.csv", "wb")
writer = csv.writer(map_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

single_females = female_population
single_males = male_population
couples = []

for generation in range(1,51):
    average_age_difference = 0.0
    pairings = 0.0
    born = 0.0
    died = 0.0

    single_females.sort(key=lambda x: x.age, reverse=False)
    single_males.sort(key=lambda x: x.age, reverse=False)

    # Pairing Phase
    for female in single_females:
        for male in single_males:
            if female.get_preference(male) == 1 and male.get_preference(female) == 1:
                couples.append([female, male])

                single_females.remove(female)
                single_males.remove(male)

                average_age_difference += abs(female.age - male.age)
                pairings += 1.0

                break

    # Reproduction Phase
    for couple in couples:
        chance = (couple[0].reproduce() + couple[1].reproduce()) / 2.0
        if chance > random.uniform(0.0, 1.0):
            born += 1

            if 0.5 > random.uniform(0.0, 1.0):
                single_females.append(Person("f", 0, couple[0].genetics))
            else:
                single_males.append(Person("m", 0, couple[1].genetics))
    # Aging single people

    for female in single_females:
        female.age_person()
        if not person.alive:
            died += 1.0
            female_population.remove(female)

    for male in single_males:
        male.age_person()
        if not person.alive:
            died += 1.0
            male_population.remove(male)

    # Aging couples and returning widows to pool
    for couple in couples:

        couple[0].age_person()
        if not couple[0].alive:
            died += 1.0

        couple[1].age_person()
        if not couple[1].alive:
            died += 1.0

        if not couple[0].alive and couple[1].alive:
            male_population.append(couple[1])
            couples.remove(couple)
        elif couple[0].alive and not couple[1].alive:
            female_population.append(couple[0])
            couples.remove(couple)
        elif not couple[0].alive and not couple[1].alive:
            couples.remove(couple)


    if pairings > 0:
        av_age_difference = average_age_difference/pairings
    else:
        av_age_difference = 0.0
    # for couple in couples:
    #     print(couple)
    total_population_size = len(single_females) + len(single_males) + (len(couples))
    print("Generation: " + str(generation).zfill(3)
    + ", Population Size: " + str(total_population_size).zfill(6)
    + ", Born: " + str(born).zfill(6)
    + ", Died: " + str(died).zfill(6)
    + ", Change: " + str(born-died).zfill(6)
    + ", Number of couples: " + str(len(couples)).zfill(4)
    + ", Average age difference: " + str(av_age_difference).zfill(6))
    writer.writerow([str(generation).zfill(3), str(total_population_size).zfill(6), str(born).zfill(6), str(died).zfill(6), str(born-died).zfill(6), str(len(couples)).zfill(4), str(av_age_difference).zfill(6)])

map_output.close()

### Begin MYP simulation ###
print("### Begin MYP Simulation ###")

myp_output = open("myp_simulation_output.csv", "wb")
writer = csv.writer(myp_output, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

single_females = female_population
single_males = male_population
couples = []

for male in single males:


for generation in range(1,51):
    average_age_difference = 0.0
    pairings = 0.0
    born = 0.0
    died = 0.0

    single_females.sort(key=lambda x: x.age, reverse=False)
    single_males.sort(key=lambda x: x.age, reverse=False)

    # Pairing Phase
    for female in single_females:
        for male in single_males:
            if female.get_preference(male) == 1 and male.get_preference(female) == 1:
                couples.append([female, male])

                single_females.remove(female)
                single_males.remove(male)

                average_age_difference += abs(female.age - male.age)
                pairings += 1.0

                break

    # Reproduction Phase
    for couple in couples:
        chance = (couple[0].reproduce() + couple[1].reproduce()) / 2.0
        if chance > random.uniform(0.0, 1.0):
            born += 1

            if 0.5 > random.uniform(0.0, 1.0):
                single_females.append(Person("f", 0, couple[0].genetics))
            else:
                single_males.append(Person("m", 0, couple[1].genetics))
    # Aging single people

    for female in single_females:
        female.age_person()
        if not person.alive:
            died += 1.0
            female_population.remove(female)

    for male in single_males:
        male.age_person()
        if not person.alive:
            died += 1.0
            male_population.remove(male)

    # Aging couples and returning widows to pool
    for couple in couples:

        couple[0].age_person()
        if not couple[0].alive:
            died += 1.0

        couple[1].age_person()
        if not couple[1].alive:
            died += 1.0

        if not couple[0].alive and couple[1].alive:
            male_population.append(couple[1])
            couples.remove(couple)
        elif couple[0].alive and not couple[1].alive:
            female_population.append(couple[0])
            couples.remove(couple)
        elif not couple[0].alive and not couple[1].alive:
            couples.remove(couple)


    if pairings > 0:
        av_age_difference = average_age_difference/pairings
    else:
        av_age_difference = 0.0
    # for couple in couples:
    #     print(couple)
    total_population_size = len(single_females) + len(single_males) + (len(couples))
    print("Generation: " + str(generation).zfill(3)
    + ", Population Size: " + str(total_population_size).zfill(6)
    + ", Born: " + str(born).zfill(6)
    + ", Died: " + str(died).zfill(6)
    + ", Change: " + str(born-died).zfill(6)
    + ", Number of couples: " + str(len(couples)).zfill(4)
    + ", Average age difference: " + str(av_age_difference).zfill(6))
    writer.writerow([str(generation).zfill(3), str(total_population_size).zfill(6), str(born).zfill(6), str(died).zfill(6), str(born-died).zfill(6), str(len(couples)).zfill(4), str(av_age_difference).zfill(6)])

map_output.close()
