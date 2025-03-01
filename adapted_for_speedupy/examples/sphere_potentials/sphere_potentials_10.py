import sys, os
from speedupy.speedupy import initialize_speedupy
from dnacc.derjaguin import calc_spheres_potential
import dnacc
from dnacc.units import nm
import numpy as np

@initialize_speedupy
def main():
    plates = dnacc.PlatesMeanField()
    L = 20 * nm
    plates.set_tether_type_prototype(sigma=0, L=L)
    ALPHA = plates.add_tether_type(plate='upper', sticky_end='alpha')
    ALPHA_P = plates.add_tether_type(plate='lower', sticky_end='alphap')
    hArr = np.linspace(1 * nm, 40 * nm, 40)
    for S in (0.1, 0.33, 0.66, 1):
        sigma = 1 / (S * L) ** 2
        plates.tether_types[ALPHA]['sigma'] = sigma
        plates.tether_types[ALPHA_P]['sigma'] = sigma
        for betaDeltaG0 in range(-7, -2):
            plates.beta_DeltaG0['alpha', 'alphap'] = betaDeltaG0
            temp3 = [plates.at(h) for h in hArr]
            betaFPlate = [h.free_energy_density for h in temp3]
            with open('plates-S%0.2f-G%.1f.dat' % (S, betaDeltaG0), 'w') as f:
                temp1 = '\t'
                f.write(temp1.join(['h / L', 'F_rep (kT/L^2)', 'F_att (kT/L^2)', 'F_plate (kT/L^2)']) + '\n')
                for (h, V) in zip(hArr, betaFPlate):
                    temp4 = plates.at(h)
                    betaFRep = temp4.rep_free_energy_density
                    betaFAtt = V - betaFRep
                    f.write('%.7g\t%.7g\t%.7g\t%.7g\n' % (h / L, betaFRep / (1 / L ** 2), betaFAtt / (1 / L ** 2), (betaFRep + betaFAtt) / (1 / L ** 2)))
            for R in (6.7, 25.0):
                betaFSphere = calc_spheres_potential(hArr, betaFPlate, R * L)
                with open('spheres-R%.1f-S%0.2f-G%.1f.dat' % (R, S, betaDeltaG0), 'w') as f:
                    temp2 = '\t'
                    f.write(temp2.join(['h / L', '[ignore: F_rep (kT)]', '[ignore: F_att (kT)]', 'F_sphere (kT)']) + '\n')
                    for (h, V) in zip(hArr, betaFSphere):
                        f.write('%.7g\t%.7g\t%.7g\t%.7g\n' % (h / L, 0, 0, V))
main()