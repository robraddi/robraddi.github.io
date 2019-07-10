#!/usr/bin/env python

import numpy as np
import itertools

class Simulation(object):
    """2D — Molecular Dynamics Simulation.

    NOTE on UNITS:
    Units are in LJ, where all quantities are unitless.
    mass, sigma, epsilon, and the Boltzmann constant = 1
    """

    def __init__(self, nSteps, dim=2, row_col=[5,10],
            mass=1., dt=1., temp=1., thermostat="Berendsen",
            print_frequency=1.,lattice="hex", lattice_spacing=[1.1,1.1],
            traj_format="lammps", write_frequency=1000,
            store_frequency=1, LJ_cutoff=2.5,
            v_initial=1., testing=False, verbose=True):

        self.verbose=bool(verbose)
        self.testing=bool(testing)
        self.print_frequency = print_frequency

        self.nSteps = int(nSteps)           #
        #self.nParticles = int(nParticles)   #
        self.epsilon = 1.                   #
        self.sigma = 1.                     #
        self.LJcutoff = LJ_cutoff*self.sigma
        self.kb = 1.                        #
        self.dt = dt                         # t (epsilon / m / sigma^2)^1/2
        self.t = 0.                         # t (epsilon / m / sigma^2)^1/2
        self.mass = mass                    # m
        self.row_col = row_col
        self.spacingx = lattice_spacing[0]   #self.spacingx = 1.120#1.214  #1.125
        self.spacingy = lattice_spacing[1]   #self.spacingy = 2.14#1.214  #0.969114

        self.dim = dim
        self.temperature=temp               # T Kb / epsilon
        self.thermostat=thermostat

        self.virialAccumulator = 0
        self.potential_energyAccumulator = 0

        self.density = self.mass/self.sigma**self.dim
        # energy (E / epsilon)

        # set initial positions and velocities
        self.q = self._set_positions(lattice)     # x / sigma
        self.v_initial = v_initial
        self.v = self._set_velocities()             # v tau / sigma

        self.forces = self.get_forces(self.q)       # f sigma / epsilon

        self.Traj = Trajectory(
                box_lengthX=self.boxlengthX,
                box_lengthY=self.boxlengthY,
                write_frequency=write_frequency,
                traj_format=traj_format)

        self.integrator = Integator(integrator="vv")

        ## Temporary fix for scaling velocities #TODO: make this user friendly
        self.scale_velocities_frequency = 10

    def LJ(self, r):
        """The Energy of Lennard-Jones Potential.

        :math:  r`V_{LJ} = 4\epsilon[(\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^{6}],`

        :param np.array r: the radial distance (units: x/sigma)
        :var float sigma: Distance at which the inter-particle potential is zero (units: None)
        :var float epsilon: the depth of the potential well (units: None)
        :return np.array: array of energies (units: E/epsilon)"""

        result = 4.*self.epsilon*( (self.sigma/r)**12. - (self.sigma/r)**6. )
        return result

    def dU(self, r):
        """The derivative of Lennard-Jones Potential.

        :param np.array r: the radial distance (units: x/sigma)
        :var float sigma: Distance at which the inter-particle potential is zero (units: None)
        :var float epsilon: the depth of the potential well (units: None)
        :return np.array: array of energies (units: E/epsilon)"""

        result = 4.*(self.epsilon/r)*(-12.*(self.sigma/r)**12 + 6.*(self.sigma/r)**6)
        return result


#    def _set_positions(self,lattice=True):
#        """Initialize random positions of N particles.
#
#        :math: r`r_{x,y} = \left( \begin{array}{c}{r_{x1} \hspace{0.25cm} r_{y1}} \\ {r_{x2} \hspace{0.25cm} r_{y2}} \\ { \vdots}\hspace{0.25cm}{ \vdots} \\ {r_{\mathrm{xN}} \hspace{0.25cm} r_{\mathrm{yN}}}\end{array}\right)`
#        :return np.array: array of positions for N particles  (units: m)"""
#
#        S = self.half_edge_length
#        S_range = np.array(range(-S,S+1),dtype=float)
#
#        if lattice=="bcc":
#            self.q = np.array(list(itertools.product(S_range, repeat=self.dim)),dtype=float)
#            self.nParticles = len(self.q)
#            self.boxlengthX = S*1.14
#            self.boxlengthY = S*1.14
#
#            #self.half_edge_length = self.boxlength/2.
#            print("Lattice Spacing in x,y = %s"%(np.abs(self.q[0] - self.q[1])))
#
#        if lattice=='hex':
#            xv,yv = np.meshgrid(S_range, S_range, sparse=False, indexing='ij')
#            X,Y = [],[]
#            for i in range(len(S_range)):
#                for j in range(len(S_range)):
#                    X.append(xv[j,i])
#                    if (j%2)!=0:
#                        yv[j,i] += 0.5
#                    Y.append(yv[j,i])
#            q = []
#            for i in range(len(Y)):
#                q.append([X[i],(Y[i] - (max(Y) + min(Y))/2.)])
#
#            self.spacingx = 1.120#1.214  #1.125
#            self.spacingy = 2.14#1.214  #0.969114
#            self.q = np.array(q)
#            self.q[:,0], self.q[:,1] = self.q[:,0]*self.spacingx, self.q[:,1]*self.spacingy
#            self.nParticles = len(self.q)
#            self.boxlengthX = self.q[-1][0]*1.2
#            self.boxlengthY = self.q[-1][-1]*1.2
#            #self.half_edge_length = self.boxlength/2.
#            print("Lattice Spacing in x,y = %s,%s"%(self.spacingx,self.spacingy))
#
#        return self.q


    def _set_positions(self,lattice=True):
        """Initialize random positions of N particles.

        :math: r`r_{x,y} = \left( \begin{array}{c}{r_{x1} \hspace{0.25cm} r_{y1}} \\ {r_{x2} \hspace{0.25cm} r_{y2}} \\ { \vdots}\hspace{0.25cm}{ \vdots} \\ {r_{\mathrm{xN}} \hspace{0.25cm} r_{\mathrm{yN}}}\end{array}\right)`
        :return np.array: array of positions for N particles  (units: m)"""

        # Get the number of rows and columns that make up the dimensions of the box
        lenX,lenY = self.row_col[0], self.row_col[1]

        # 1-D box if lengeth of x or y eq 1
        if lenX == 1:
            X_range = np.array([1],dtype=float)
            particle_wall_spacingX = 1.0
        # Anything else is 2-D
        else:
            if lenX%2: # if x is odd
                X_range = np.array(range(-int(lenX/2),int(lenX/2)+1),dtype=float)
            else:      #if x is even
                X_range = np.array(range(-int(lenX/2),int(lenX/2)),dtype=float)

        if lenY == 1:
            Y_range = np.array([1],dtype=float)
            particle_wall_spacingY = 1.0
        else:
            if lenY%2:
                Y_range = np.array(range(-int(lenY/2),int(lenY/2)+1),dtype=float)
            else:
                Y_range = np.array(range(-int(lenY/2),int(lenY/2)),dtype=float)

        xv,yv = np.meshgrid(X_range, Y_range, sparse=False, indexing='ij')

        if lattice=="bcc":
            q = []
            for i in range(len(Y_range)):
                for j in range(len(X_range)):
                    q.append([xv[j,i],yv[j,i]])
            extra = 0.5

        if lattice=='hex':
            X,Y = [],[]
            for i in range(len(Y_range)):
                for j in range(len(X_range)):
                    X.append(xv[j,i])
                    if (j%2)!=0:
                        yv[j,i] += 0.5
                    Y.append(yv[j,i])
            q = []
            for i in range(len(Y)):
                q.append([X[i],(Y[i] - (max(Y) + min(Y))/2.)])
            extra = 0

        self.q = np.array(q)

        # Map the coordinates to the origin by shift
        self.q[:,0] += abs(self.q[0,0])
        self.q[:,1] += abs(self.q[0,1])
        self.nParticles = lenX*lenY

        # Scale the spacing in between particles
        self.q[:,0], self.q[:,1] = self.q[:,0]*self.spacingx, self.q[:,1]*self.spacingy

        self.boxlengthX = abs(self.q[-1,0])+0.4+extra
        self.boxlengthY = abs(self.q[-1,1])+0.8+extra

        print("# of Particles = %s"%self.nParticles)
        print("BOX DIMENSIONS = (%s,%s)"%(self.boxlengthX,self.boxlengthY))
        print("Lattice Spacing in x,y = %s,%s"%(self.spacingx,self.spacingy))
        #print(self.q)
        #exit(1)
        return self.q



    def _set_velocities(self):
        """Initialize random velocities of N particles.

        :math: r`v_{x,y} = \left( \begin{array}{c}{v_{x1} \hspace{0.25cm} v_{y1}} \\ {v_{x2} \hspace{0.25cm} v_{y2}} \\ { \vdots}\hspace{0.25cm}{ \vdots} \\ {v_{\mathrm{xN}} \hspace{0.25cm} v_{\mathrm{yN}}}\end{array}\right)`
        :return np.array: array of velocities for N particles  (units: `m s^-1`)"""

        # Create Nx2 array of random velocities for N particles
        self.v = self.v_initial*np.transpose(np.array(
            [np.random.normal(size=self.nParticles),
            np.random.normal(size=self.nParticles)]))
        vsum = self.v.sum(axis=0)/self.nParticles

        self.v[:,0] -= vsum[0]
        self.v[:,1] -= vsum[1]

        # TODO: make sure that velocities go to the Boltzmann distribution
        # vi = np.sqrt(mi/2pikT)np.exp(-mv**2/(2kT))

        # Scale the velocities according to the temperature
        self.scale_velocities()
        return self.v

    def scale_velocities(self):
        """Return velocities for N particles.

        :math: r`\sqrt{\frac{k_{b}T}{m}}`

        Isokinetics thermostat:
        :math: r`v(t) = v(t)\sqrt{\frac{SetTemp}{T}}`

        Berendsen thermostat:
        :math: r`v_{x,y} = v_{x,y}*\left[1+\frac{\nu \Delta t}{\tau_{T}}\left\{\frac{T_{0}}{T\left(t-\frac{1}{2} \Delta t\right)}-1\right\}\right]^{1 / 2}`

        :return np.array: (units: `m s^-1`)"""

        #pass
        # measure the instantaneous temperature
        T_now = self.get_temperature()

        if self.thermostat == "Isokinetics":
            # Scale the velocities to achieve the desired temperature
            self.v = self.v*np.sqrt(self.temperature/T_now)

        if self.thermostat == "Berendsen":
            # (the larger the rise_time #, the longer to achieve desired temp)
            rise_time = 0.01 # 0.01, 0.1, 1.0, 10
            self.v = self.v*np.sqrt(1 + (self.dt/rise_time)*(self.temperature/T_now - 1.))

        return self.v



    def get_temperature(self):
        """Return the temperature for N particles each with 3 degrees of
        freedom (DOF).

        :math: r`E_{\mathrm{kinetic}} = \left\langle\frac{1}{2} m v_{\alpha}^{2}\right\rangle=\frac{1}{2} k_{B} T`
        :math: r`T=\frac{2 k_{\mathrm{B}} E_{\mathrm{kinetic}}}{3 N}`
        :return float: temperature (units: K)"""

        #self.T = 0.5*self.get_kinetic_energy()/(self.nParticles*self.kb)
        self.T = self.dim*0.5*self.get_kinetic_energy()/(self.nParticles*self.kb)
        return self.T

    def get_pressure(self):
        """Returns the pressure for the system."""

        #self.Pressure = self.density*self.kb*self.T+1/(self.dim)*sum([self.forces[i,0]+self.forces[i,1] for i in range(len(self.forces))])/self.nParticles
        step = self.step
        if step == 0:
            step = 0.1
        Virial = self.virialAccumulator/step
        vol = self.boxlengthX*self.boxlengthY*abs(-1e-1 - 1e-1)
        self.Pressure = self.nParticles*self.kb*self.T/vol+1./(self.dim*vol)*Virial/(self.nParticles*self.T)
        return self.Pressure


    def check_pbc(self,q):
        """ Check the periodic boundary conditions for all the particles"""

        for i in range(self.nParticles):
            # For x
            if (q[i,0] >= self.boxlengthX):
                q[i,0] -= self.boxlengthX
            elif (q[i,0] <= 0.0):
                q[i,0] += self.boxlengthX
            # For y
            if (q[i,1] >= self.boxlengthY):
                q[i,1] -= self.boxlengthY
            elif (q[i,1] <= 0.0):
                q[i,1] += self.boxlengthY
        return q


    def get_kinetic_energy(self):
        """Return the kinetic energy of N particles

        :math: r`E_{k i n}=\frac{1}{2} \sum_{i=1}^{N} m_{i} v_{i}^{2}`
        """

        result = 0.5*(self.mass*self.v**2.).sum()
        return result


    def get_avg_kinetic_energy(self):
        """Return the average kinetic energy of the system."""

        return self.dim*self.nParticles*self.kb*self.T/2.


    def get_forces(self,q):
        """Returns the forces for N particles.

        :math: r`r_{ij} = \sum_{i=1}^{N}\sum_{j=1\\j\neq i}^{N} |r_{i} - r_{j}|`
        :math: r`F = - (\frac{dU(r)}{dr})`
        :math: r`F_{x} = - (\frac{x}{r})(\frac{dU(r)}{dr})\hspace{0.25cm};\hspace{0.25cm} F_{y} = - (\frac{y}{r})(\frac{dU(r)}{dr})`
        :math: r` F = F_{x} + F_{y}`

        :return np.array: the forces for N particles (units: N)"""

        # Reset the forces
        self.forces = np.transpose([np.zeros(self.nParticles),
                    np.zeros(self.nParticles)])

        self.potential_E = 0
        #self.virialAccumulator = 0
        for i in range(0,self.nParticles-1):
            for j in range(i+1,self.nParticles):
                # x and y components of the ith and jth particle
                RX, RY = q[j,0]-q[i,0], q[j,1]-q[i,1]

                #print(RX,RY)
                #exit(1)
                # Conditional statements for PBC
                if (RX > self.boxlengthX/2.):
                    RX -= self.boxlengthX
                elif (RX < -self.boxlengthX/2.):
                    RX += self.boxlengthX
                if (RY > self.boxlengthY/2.):
                    RY -= self.boxlengthY
                elif (RY < -self.boxlengthY/2.):
                    RY += self.boxlengthY

                R2 = np.sqrt(RX*RX + RY*RY)
                # Implementing the Lennard-Jones cutoff and subtracting
                # the associated energy
                if R2 <= self.LJcutoff:
                    #F = self.dU(R2)
                    U = self.LJ(R2) - self.LJ(self.LJcutoff)
                else:
                    #F = 0
                    U = 0

                F = self.dU(R2)
                if F > 20000:
                    print(i,j,R2,F)
                    print("Simulation Exploded!...")
                    exit(1)

                self.forces[i,0] = self.forces[i,0]+F*RX/R2
                self.forces[i,1] = self.forces[i,1]+F*RY/R2
                self.forces[j,0] = self.forces[j,0]-F*RX/R2
                self.forces[j,1] = self.forces[j,1]-F*RY/R2
                self.potential_E += U
                self.virialAccumulator += F*RX/R2+F*RY/R2


        self.T = self.get_temperature()

        #if self.Pressure < 1e1:
        #    scaling_factor = 1e1
        #    self.forces = self.forces*scaling_factor

        #sigma = 3.4e-10    #3.4 × 10-10 m
        #epsilon = 1.65e-21 #1.65 × 10-21 J
        #print(self.forces*epsilon/sigma)
        #print(self.forces)

        #if max(np.transpose(self.forces)[0]) > 1e2:
        #    print(self.forces)
        #    exit(1)
        return self.forces




    def integrate(self):
        """Integration of Newton's laws.

        :param r: the positions of N particles (units: m)
        :param v: the velocities of N particles (units: m s^-1)
        :returns: r, v, T, total_energy    (units:   )
        """

        print("NUMBER OF PARTICLES = %s"%self.nParticles)
        if self.verbose:
            print("""
Step #\tTime,t \t\tTemperature,T \t\tPair Potential, U\tTotal Energy, E\t\tPressure, P""")


        # Keep track of the number of trajectories written
        ntraj = 0
        for step in range(self.nSteps):
            self.step = step



            # Get the acceleration from the forces
            a = self.forces/self.mass

            if step%self.scale_velocities_frequency == 0:
                if self.testing:
                    print("Scaling velocities...")
                self.scale_velocities()

            # Velocity Verlet algorithm:  ri(t+dt) = ri(t) + vi(t)dt + (a(t)/2) dt^2
            qnew = self.q + self.v*self.dt + 0.5*a*self.dt**2

            # New velocities (requires calculating the accelerations at the new positions)
            Fnew = self.get_forces(qnew)
            anew = Fnew/self.mass
            vnew = self.v + 0.5*(a + anew)*self.dt

            if self.testing:
                for i in range(len(vnew)):
                    print("qnew = %s m"%qnew[i])
                for i in range(len(vnew)):
                    print("vnew = %s m s^-1"%vnew[i])

            self.KE = self.get_kinetic_energy()

            # Make sure that we get the initial step
            if step == 0:
                self.Traj.trajectory.append(np.array([self.q, self.v]))

            # Update positions and velocities after correcting Kinetic Energy
            self.q, self.v, self.forces = qnew, vnew, Fnew

            # periodic boundary conditions
            self.q = self.check_pbc(self.q)

            if self.testing:
                print("qnew = %s m "%self.q)
                print("vnew = %s m s^-1"%self.v)

            # Update the temperature
            self.get_temperature()

            self.get_pressure()

            # Update the total energy
            #self.total_energy = (self.potential_E+self.KE)

            #TODO: Make this nicer. Way too ugly
            self.potential_energyAccumulator += self.potential_E
            STEP = self.step
            if self.step == 0:
                STEP = 0.1
            self.potential_E = self.potential_energyAccumulator/STEP


            self.total_energy = (self.potential_E+self.get_avg_kinetic_energy())

            # Energy per particle for N particles
            self.energy_per_particle = self.total_energy/self.nParticles

            if self.testing:
                print('Total Energy = %s kJ'%self.total_energy)

            if self.verbose:
                if (step)%self.print_frequency == 0:
                        #print("{: >20} {: >20} {: >20}".format(*row))
                    print("%i\t%0.4e\t%0.4e\t\t%0.4e\t\t%0.4e\t\t%0.4e\t"%(
                        step, self.t, self.T, self.potential_E,
                        self.total_energy,self.Pressure))
            #exit(1)
            self.t += self.dt

            if step%self.Traj.store_frequency == 0:
                self.Traj.trajectory.append(np.array([self.q, self.v]))

            if (step+1) == self.Traj.write_frequency:
                print("\nWriting Trajectory: traj_%s.%s\n"%(
                    ntraj,self.Traj.traj_format))
                if self.Traj.traj_format=="npz":
                    self.Traj.write_npy(outfilename="traj_%s.npz"%(ntraj))
                if self.Traj.traj_format=="lammps":
                    self.Traj.write_lammps("traj_%s.lammpstraj"%(ntraj),
                            nParticles=self.nParticles)
                if self.Traj.traj_format=="xyz":
                    self.Traj.write_xyz("traj_%s.xyz"%(ntraj),
                            nParticles=self.nParticles)
                ntraj += 1

        #return self.q, self.v, self.T, self.total_energy

class Integator(object):
    """Class of integrators to use for the simulation"""


    def __init__(self, integrator="vv"):
        """Types of integrators:
            vv = Velocity Verlet
            brown = Over-damped Langevin (Brownian Dynamics) """

        self.integrator = integrator


    def velocity_verlet(self):
        """The velocity Verlet integrator.

        :math: r`r(t+\Delta t)=r(t)+v(t) \Delta t+\frac{F_{LJ}(t)}{2 m} \Delta t^{2}`
        :math: r`r_{i}(t+\Delta t) = r_{i}(t) + v_{i}(t) + (a(t)/2) \Delta t^2 `
        :math: r`a(t+\Delta t) = F_{LJ}(r(t+\Delta t))/m`
        :math: r`v(t+\Delta t) = v(t) + 0.5 (a(t+\Delta t)+a(t))\Delta t `

        """

        pass

    def brownian_dynamics(self):
        """Coupled over-damped Langevin. Here, the equations of motion are
        non-dimensionalized using time as `\tau = \sigma^{2}/D`, where
        `\sigma` and `k_{B}T` as basic units of length and energy.

        :math: r`\dot{\boldsymbol{r}}_{i}=\frac{1}{\gamma} \boldsymbol{F}_{\mathrm{LJ}}\left(\left\{\boldsymbol{r}_{i}\right\}\right)+v_{\mathrm{p}} \hat{\boldsymbol{v}}_{i}+\sqrt{2 D} \boldsymbol{\eta}_{i}^{\mathrm{T}}`
        :math: r`\hat{v}_{i}=\left(\cos \theta_{i}, \sin \theta_{i}\right)`
        :math: r`\hat{\theta}=\sqrt{2 \mathcal{D}_{r}} \eta_{i}^{R}`
        :math: r`D_{r} = 3D/\sigma^{2} \text{ at a low Reynolds number.}`
        :math: r`D=\frac{k_{\mathrm{B}} T}{\gamma}`

        """
        pass



class Trajectory(object):
    """Container class to store the simulation trajectory data. Writes trajectory
    to various file types."""

    def __init__(self, box_lengthX, box_lengthY, store_frequency=1, write_frequency=1000,
            traj_format="lammps"):
        self.trajectory = []
        self.store_frequency = store_frequency
        self.write_frequency = write_frequency
        self.traj_format = traj_format
        self.box_lengthX = box_lengthX
        self.box_lengthY = box_lengthY


    def write_npz(self, outfilename):
        """Writes a compact file of several arrays into binary format.
        Standardized: Yes ; Binary: Yes; Human Readable: No;

        :param str outfilename: name of the output file
        :return: numpy compressed filetype
        """

        np.savez_compressed(outfilename, self.trajectory)

    def read_npz(self,filename):
        """Reads a numpy compressed filetype(*.npz) file"""

        loaded = np.load(filename)
        files = loaded.files
        for step in loaded["arr_0"]:
            print(step)


    def write_xyz(self, filename, nParticles, atom="Ar"):
        """Writes a XYZ trajectory (*.xyz) file"""

        import os
        new_dir = "XYZ"
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        curDir = os.path.abspath(os.getcwd())
        t = 0
        first, last = "first", "last"

        for step in self.trajectory:
            #Natoms = '%s\n'%(nParticles)
            traj_name = filename.split("_")[0]+str("_%s"%t)+".xyz"
            fullpath = curDir+"/"+new_dir+"/"+traj_name
            if t==0:
                first = fullpath
            #print(fullpath)
            with open(fullpath, 'w') as file:
                file.seek(0)

#                comment = '''
#Lattice="%s %s 0.0 %s %s 0.0 %s %s 0.0" Properties=species:S:1:pos:R:3:select:I:1 Time=%s\n'''%(
#        -float(self.box_length), float(self.box_length),
#        -float(self.box_length), float(self.box_length),
#        -1e-1, 1e-1,
#        t)
                #file.writelines(Natoms)
                #file.writelines("liquid argon\n")
                #part = 1
                for i in range(len(step[0])):
                    line = "%s\t%s\t%s\t%s\n"%(atom, step[0][i][0], step[0][i][1], 0.0)
                    file.writelines(line)
                #    part += 1
                t += 1
                file.writelines("\n")
                file.truncate()

        last = fullpath





    def write_lammps(self, filename, nParticles, atom="Ar"):
        """Writes a XYZ trajectory (*.xyz) file"""

        #sigma = 3.4e-10
        traj_name = filename.split(".")[0]+".lammpstrj"
        with open(traj_name, 'w') as file:
            file.seek(0)
            t = 0
            for step in self.trajectory:
                file.writelines("ITEM: TIMESTEP\n%s\n"%t)
                file.writelines("ITEM: NUMBER OF ATOMS\n%s\n"%nParticles)
                file.writelines("ITEM: BOX BOUNDS pp pp pp\n")
                file.writelines("%s %s\n%s %s\n%s %s\n"%(
                    0.0, float(self.box_lengthX),
                    0.0, float(self.box_lengthY),
                    -1e-1, 1e-1))
                file.writelines("ITEM: ATOMS id type x y z vx vy vz \n")
                part = 1
                for i in range(len(step[0])):
                    line = "%s %s %s %s %s %s %s %s\n"%(
                            part, 1, step[0][i][0], step[0][i][1], 0,
                            step[1][i][1], step[1][i][1], 0)
                    file.writelines(line)
                    part += 1
                t += 10
            file.truncate()


class Analysis(object):
#class Analysis(Trajectory):    ????
    """Analyze the trajectory"""

    def __init__(self):
        pass
        self.self = 1


    def get_RDF(self):
        """Returns the radial distribution function, g(r)

        :math: r`g(r) = \rho_{}`

        """

        pass

        for i in range(len(self.nParticles)):
            dx -= lx*nint(dx/lx)
            dy -= ly*nint(dy/ly)
            dz -= lz*nint(dz/lz)
            d = math.sqrt(dx*dx + dy*dy + dz*dz)
            bin = int(d/dr)
            if (bin < len(counts)):
                counts[bin] += 1


    def get_nearest_neighbors(self):
        """Returns the nearest neighbors, n(r) from integrating g(r)

        :math: r`n(r) = 4\pi \rho \int_{0}^{r}r^{2}g(r)dr`

        """

        pass

        gr = self.get_RDF()

        intgr = 0.0
        for i in range(0,len(counts)):
            r = dr*(i+0.5)
            f = 4. * np.pi * dens2 * r*r * dr
            u = float(counts[i])/(nstructs*n1)
            intgr += u
            print(r, u/f, intgr)

    def get_PMF(self):
        """Returns the potential mean force (PMF), w(r).

        :math: r`w(r) = -k_{B}T ln(g(r))`

        """

        pass
        gr = self.get_RDF()
        result = -self.kb*self.T*np.log(gr)
        return result



if __name__ == "__main__":

    ###########################################################################
    #################  MAIN  ##################################################
    ###########################################################################
    ########## Initialize Simulation ##########################################
    NSTEPS = 1000
    #row_col = [12,4]
    row_col = [24,8]
    lattice_spacing=[0.7,1.75]  # if hex
    #lattice_spacing=[1.42,1.42]   # if bcc
    MD = Simulation(nSteps=NSTEPS, row_col=row_col, lattice="hex",#"bcc""hex",
            lattice_spacing=lattice_spacing,
            print_frequency=100, dt=0.005, v_initial=2.5, temp=1.65,
            write_frequency=NSTEPS, store_frequency=1,traj_format="xyz")#"lammps")
    ########## Integrate Equations of Motion ##################################
    # Verlet Algorithm
    MD.integrate()










