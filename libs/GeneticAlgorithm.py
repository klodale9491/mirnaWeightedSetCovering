# Coder : Alessio Giorgianni

from random import randint
from random import random
import matplotlib.pyplot as plt
import csv

class GeneticAlgorithm:

    def __init__(self,
                 fitness_fnc,     # Funzione di fitness custom
                 fitness_cmp_fnc, # Permette di confrontare 2 valori di fitness. Torna 1 se la prima fit
                                  # e' migliore della seconda, 0 altrimenti.
                 min_val_gene,    # Min valore per del gene (nel caso binario, 0)
                 max_val_gene,    # Max valore per del gene (nel caso binario, 1)
                 pop_len,         # Cardinalita della popolazione iniziale
                 sol_len         # Numero di geni(bit) del cromosoma(soluzione)
                 ):
        self.fitness_fnc = fitness_fnc
        self.fitness_cmp_fnc = fitness_cmp_fnc
        self.min_val_gene = min_val_gene
        self.max_val_gene = max_val_gene
        self.pop_len = pop_len
        self.sol_len = sol_len
        self.surv_len = int(pop_len / 2)
        # Popolazione delle soluzioni
        self.population = []
        self.fitness = []
        self.random_population()
        # Performarce algoritmo
        self.times = []
        self.scores  = [] # np coperte
        self.scores2 = [] # nr coperte
        self.scores3 = [] # grandezza subset


    # Plotta le performance dell'algoritmo genetico
    def plot_performance(self):
        # Plotto i grafici
        fig = plt.figure()

        plt.subplot(2, 2, 1)
        plt.xlabel('generations')
        plt.ylabel('np_covered')
        plt.plot(self.times, self.scores)
        plt.grid(True)

        plt.subplot(2, 2, 2)
        plt.xlabel('generations')
        plt.ylabel('nr_covered')
        plt.plot(self.times, self.scores2)
        plt.grid(True)

        plt.subplot(2, 2, 3)
        plt.xlabel('generations')
        plt.ylabel('subsets')
        plt.plot(self.times, self.scores3)
        plt.grid(True)

        plt.show()


    # Imposta il dataset ed altre strutture dati utili
    def set_params(self, par1, par2):
        self.mirna_data = par1
        self.up_c = par2


    # Calcola il numero di proteine up regolate della soluzione
    def get_np_covered(self, solution):
        # Proteine up_regulated che ho coperto
        prots_up = [0 for i in range(self.up_c)]
        for i in range(len(solution)):
            if solution[i] == 1:
                for x in range(len(prots_up)):
                    prots_up[x] |= self.mirna_data[i][x]
        return len(prots_up) - prots_up.count(0)


    # Calcola il numero di proteine normali della soluzione
    def get_nr_covered(self, solution):
        # Proteine up_regulated che ho coperto
        prots_nr = [0 for i in range(self.up_c, len(self.mirna_data[0]))]
        for i in range(len(solution)):
            if solution[i] == 1:
                for y in range(len(prots_nr)):
                    prots_nr[y] |= self.mirna_data[i][y]
        return len(prots_nr) - prots_nr.count(0)


    def save_data(self, run_id, path):
        with open(path + "run_" + str(run_id) + ".csv", "wb") as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for i in range(len(self.scores)):
                spamwriter.writerow([self.scores[i], self.scores2[i], self.scores3[i]])
            csvfile.close()


    # Trova la migliore fitness tornando la soluzione
    # ed il suo valore
    def fitness_bst(self,fits):
        best_index = 0
        for i in range(1, len(fits)):
            if self.fitness_cmp_fnc(fits[i], fits[best_index]) == 1:
                best_index = i
        return fits[best_index], best_index


    def fitness_bst2(self, fits):
        max_val = max(fits)
        max_idx = fits.index(max_val)
        return max_idx


    # Funzione di selezione naturale in base al metodo a TORNEO
    # Permette di far sopravvivere una percentuale delle solzioni
    def natural_selection(self, N, fit):
        new_pop = []
        for i in range(N):
            # Scelgo a caso due elementi della popolazione e verifico chi sia il migliore aggiunggendoli
            # alla nuova popolazione.
            idx_1 = randint(0, len(self.population) - 1)
            idx_2 = randint(0, len(self.population) - 1)
            if self.fitness_cmp_fnc(fit[idx_1], fit[idx_2]) == 1:
                new_pop.append(self.population[idx_1])
            else:
                new_pop.append(self.population[idx_2])
        return new_pop


    # Selezione con fitness proporzionale
    def natural_selection_2(self,N):
        new_pop = []
        fit_vec = [0 for i in range(self.pop_len + 1)]
        for i in range(1,len(fit_vec)):
            fit_vec[i] = fit_vec[i-1] + self.fitness_fnc(self.population[i-1])
        for i in range(N):
            # scelgo un numero fra 0 e fit_vec[-1](max)
            r = random() * fit_vec[-1]
            for j in range(1,len(fit_vec)):
                if r >= fit_vec[j-1] and r <= fit_vec[j] :# attenzione al minore uguale
                    id = j - 1
                    new_pop.append(self.population[id])
                    break
        return new_pop


    # Funzione che crea una popolazione random di soluzioni
    def random_population(self):
        # Creo soluzioni di taglia sol_len
        for i in range(self.pop_len):
            sol = []
            for j in range(self.sol_len):
                sol.append(randint(self.min_val_gene, self.max_val_gene))
            self.population.append(sol)


    # Funzione che effettua crossover a due punti della soluzione
    def crossover2(self, sol_1, sol_2):
        pnt_1 = randint(0,len(sol_1)-1)
        pnt_2 = randint(pnt_1, len(sol_1)-1)
        # Scambio dei geni ...
        s1 = [sol_1[0:pnt_1], sol_2[pnt_1:pnt_2], sol_1[pnt_2:]]
        s2 = [sol_2[0:pnt_1], sol_1[pnt_1:pnt_2], sol_2[pnt_2:]]
        # Ricompongo le soluzioni
        sol_1 = []
        sol_2 = []
        for i in range(3):
            for j in range(len(s1[i])):
                sol_1.append(s1[i][j])
                sol_2.append(s2[i][j])
        return sol_1 if (self.fitness_cmp_fnc(self.fitness_fnc(sol_1), self.fitness_fnc(sol_2)) == 1) else sol_2


    # Crossover ad un punto
    def crossover(self, sol_1, sol_2):
        frg_len = randint(1,len(sol_1)-1)
        s1 = [sol_1[0:frg_len], sol_2[frg_len:]]
        s2 = [sol_2[0:frg_len], sol_1[frg_len:]]
        # Ricompongo le soluzioni
        sol_1 = []
        sol_2 = []
        for i in range(2):
            for j in range(len(s1[i])):
                sol_1.append(s1[i][j])
                sol_2.append(s2[i][j])
        return sol_1 if(self.fitness_cmp_fnc(self.fitness_fnc(sol_1), self.fitness_fnc(sol_2)) == 1) else sol_2


    # Funzione che simula una mutazione puntuale di un gene
    def mutation(self, sol, pm):
        for i in range(len(sol)):
            if random() <= pm:
                sol[i] = randint(self.min_val_gene, self.max_val_gene)
        return sol


    # Main routine dell'algoritmo genetico
    def run(self, gen, pc=0.5, pm=0.05, debug=False):
        good_sol = self.population[0]
        stall = 0
        MAX_STALL = 200
        for g in range(gen):
            fitness = [self.fitness_fnc(self.population[i]) for i in range(len(self.population))]
            if debug:
                best_vec = self.fitness_bst2(fitness)
                if self.fitness_cmp_fnc(self.population[best_vec], good_sol) == 1:
                    good_sol = self.population[best_vec]
                fit_avg = sum(fitness) / len(fitness)
                best_sol = self.population[best_vec]
                np_covered = self.get_np_covered(best_sol)
                nr_covered = self.get_nr_covered(best_sol)
                subsets = sum(best_sol)
                self.times.append(g)
                self.scores.append(np_covered)
                self.scores2.append(nr_covered)
                self.scores3.append(subsets)
                '''
                if g >= 1:
                    # Verifica di stallo
                    if self.scores[g] == self.scores[g-1] and self.scores2[g] == self.scores2[g-1] and self.scores3[g] == self.scores3[g-1]:
                        stall += 1
                    else:
                        stall = 0
                if stall >= MAX_STALL:
                    return good_sol
                '''
                print('generation #', g, 'good_sol_val: ', self.fitness_fnc(good_sol), 'fit_avg: ', fit_avg, 'best np_covered: ', np_covered, 'best nr_covered: ', nr_covered, 'subsets: ', subsets)
            # Selezione naturale
            new_population = self.natural_selection(self.surv_len, fitness)
            # Crossover
            for i in range(self.surv_len):
                if pc >= random():
                    sol = self.crossover2(new_population[i], new_population[randint(0, len(new_population) - 1)])
                    new_population.append(sol)
            # Mutazioni
            for k in range(len(new_population)):
                new_population[k] = self.mutation(new_population[k], pm)
            # Aggiornamento della popolazione
            self.population[:] = new_population[:]
        return good_sol
