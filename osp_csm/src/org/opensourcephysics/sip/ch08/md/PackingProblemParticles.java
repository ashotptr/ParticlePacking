/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see:
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.sip.ch08.md;
import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;
import java.util.Scanner;

import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;

/**
 * LJParticlesApp evolves a two-dimensional system of interacting particles
 * via the Lennard-Jones potential using a Verlet ODESolver.
 *
 * getHeatCapacity method corrected based on bug report by Mike Cooke.
 *
 * @author Jan Tobochnik, Wolfgang Christian, Harvey Gould
 * @version 1.1 revised 01/14/06
 */
public class PackingProblemParticles implements Drawable, ODE {
  public double state[];
  public double ax[], ay[];
  public int N, nx, ny; // number of particles, number per row, number per column
  public double Lx, Ly;
  public double rho = N/(Lx*Ly);
  public double initialKineticEnergy;
  public int steps = 0;
  public double dt = 0.01;
  public double t;
  public double totalPotentialEnergyAccumulator;
  public double totalKineticEnergyAccumulator, totalKineticEnergySquaredAccumulator;
  public double virialAccumulator;
  public String initialConfiguration;
  public double radius = 0.5;
  Verlet odeSolver = new Verlet(this);
  public double g = 1.0;
  public double k = 1.0;
  public double damping = 0.1;

  public int n = 14;           
  public int m = 4;            

  public double epsilon = 1.0;          
  public double sigma = 1.0;           
  public double mass = 1.0;        
  public double rcut = 2.5 * sigma; 
  public double rcut2 = rcut * rcut;
  public double minX, maxX, minY, maxY;
  public double width, height, area;
  public double[] size;

  private void computeBoundingBox() {
    if (N == 0) {
      minX = maxX = minY = maxY = width = height = area = 0;
      return;
    }
    
    double x = state[0];
    double y = state[2];
    minX = x - size[0];
    maxX = x + size[0];
    minY = y - size[0];
    maxY = y + size[0];

    for (int i = 1; i < N; i++) {
      x = state[4*i];
      y = state[4*i + 2];
      double left = x - size[i];
      double right = x + size[i];
      double bottom = y - size[i];
      double top = y + size[i];

      if (left < minX) minX = left;
      if (right > maxX) maxX = right;
      if (bottom < minY) minY = bottom;
      if (top > maxY) maxY = top;
    }

    width = maxX - minX;
    height = maxY - minY;
    area = width * height;
  }

  public void readFromFile(String filename) throws IOException {
    try (Scanner sc = new Scanner(new File(filename))) {
      N = sc.nextInt();
      state = new double[1 + 4 * N];
      size = new double[N];
      for (int i = 0; i < N; i++) {
        state[4*i] = sc.nextDouble();
        state[4*i+2] = sc.nextDouble();
        size[i] = sc.nextDouble();
      }
    }
  }
  
  public void initialize() {
    N = nx*ny;
    size = new double[N];
    java.util.Random rand = initialConfiguration.equals("fixed") ? 
          new java.util.Random(12344) : new java.util.Random();

    for (int i = 0; i < N; i++) {
      size[i] = 0.3 + 0.7 * rand.nextDouble();
    }
    
    t = 0;
    rho = N/(Lx*Ly);
    resetAverages();
    state = new double[1+4*N];
    ax = new double[N];
    ay = new double[N];
    if (initialConfiguration.equals("file")) {
      try {
        readFromFile("input.txt");
      } catch (IOException e) {
        e.printStackTrace();
      }
    } else if(initialConfiguration.equals("triangular")) {
      setTriangularLattice();
    } else if(initialConfiguration.equals("rectangular")) {
      setRectangularLattice();
    } else {
      setRandomPositions();
    }
    setVelocities();
    computeAcceleration();
    odeSolver.setStepSize(dt);
  }

  public void setRandomPositions() {
    java.util.Random rand = initialConfiguration.equals("fixed") ? 
          new java.util.Random(12344) : new java.util.Random();
      
    double rMinimumSquared = Math.pow(2.0, 1.0/3.0);
    boolean overlap;
    for(int i = 0;i<N;++i) {
      do {
        overlap = false;
        state[4*i] = Lx * rand.nextDouble();
        state[4*i+2] = Ly * rand.nextDouble();
        int j = 0;
        while((j<i)&&!overlap) {
          double dx = state[4*i]   - state[4*j];
          double dy = state[4*i+2] - state[4*j+2];
          double minDist = size[i] + size[j];
          if(dx*dx+dy*dy<minDist * minDist) {
            overlap = true;
          }
          j++;
        }
      } while(overlap);
    }
  }

  public void setRectangularLattice() { 
    double dx = Lx/nx;
    double dy = Ly/ny;
    for(int ix = 0;ix<nx;++ix) { 
      for(int iy = 0;iy<ny;++iy) {
        int i = ix+iy*ny;
        state[4*i] = dx*(ix+0.5);
        state[4*i+2] = dy*(iy+0.5);
      }
    }
  }

  public void setTriangularLattice() {
    double dx = Lx/nx;
    double dy = Ly/ny;
    for(int ix = 0;ix<nx;++ix) {
      for(int iy = 0;iy<ny;++iy) {
        int i = ix+iy*ny;
        state[4*i+2] = dy*(iy+0.5);
        if(iy%2==0) {
          state[4*i] = dx*(ix+0.25);
        } else {
          state[4*i] = dx*(ix+0.75);
        }
      }
    }
  }

  public void setVelocities() {
    // double vxSum = 0.0;
    // double vySum = 0.0;
    // for(int i = 0;i<N;++i) {            // assign random initial velocities
    //   state[4*i+1] = Math.random()-0.5; // vx
    //   state[4*i+3] = Math.random()-0.5; // vy
    //   vxSum += state[4*i+1];
    //   vySum += state[4*i+3];
    // }
    // // zero center of mass momentum
    // double vxcm = vxSum/N; // center of mass momentum (velocity)
    // double vycm = vySum/N;
    // for(int i = 0;i<N;++i) {
    //   state[4*i+1] -= vxcm;
    //   state[4*i+3] -= vycm;
    // }
    // double v2sum = 0; // rescale velocities to obtain desired initial kinetic energy
    // for(int i = 0;i<N;++i) {
    //   v2sum += state[4*i+1]*state[4*i+1]+state[4*i+3]*state[4*i+3];
    // }
    // double kineticEnergyPerParticle = 0.5*v2sum/N;
    // double rescale = Math.sqrt(initialKineticEnergy/kineticEnergyPerParticle);
    // for(int i = 0;i<N;++i) {
    //   state[4*i+1] *= rescale;
    //   state[4*i+3] *= rescale;
    // }
    for (int i = 0; i < N; i++) {
      state[4*i+1] = 0;
      state[4*i+3] = 0;
    }
  }

  public double getMeanTemperature() {
    return totalKineticEnergyAccumulator/(N*steps);
  }

  public double getMeanEnergy() {
    return totalKineticEnergyAccumulator/steps+totalPotentialEnergyAccumulator/steps;
  }

  public double getMeanPressure() {
    double meanVirial;
    meanVirial = virialAccumulator/steps;
    return 1.0+0.5*meanVirial/(N*getMeanTemperature());
  }

  public double getHeatCapacity() {
     double meanTemperature = getMeanTemperature();
     double meanKineticEnergySquared = totalKineticEnergySquaredAccumulator/steps;
     double meanKineticEnergy = totalKineticEnergyAccumulator/steps;
     
     double sigma2 = meanKineticEnergySquared-meanKineticEnergy*meanKineticEnergy;
     double denom = 1.0 - sigma2/(N*meanTemperature*meanTemperature);
     return N/denom;
  }

  public void resetAverages() {
    steps = 0;
    virialAccumulator = 0;
    totalPotentialEnergyAccumulator = 0;
    totalKineticEnergyAccumulator = 0;
    totalKineticEnergySquaredAccumulator = 0;
  }
  
  public void computeAcceleration() {
    for (int i = 0; i < N; i++) {
      ax[i] = 0;
      ay[i] = 0;
    }
    totalPotentialEnergyAccumulator = 0;
    virialAccumulator = 0;

    for (int i = 0; i < N - 1; i++) {
      for (int j = i + 1; j < N; j++) {
        double dx = state[4 * i] - state[4 * j];
        double dy = state[4 * i + 2] - state[4 * j + 2];
        double r2 = dx * dx + dy * dy;
        if (r2 == 0) continue; 

        double r = Math.sqrt(r2);
        double sigma_ij = size[i] + size[j];
        double f_repulsive = epsilon * n * Math.pow(sigma_ij, n) / Math.pow(r, n + 1);
        double f_attraction = epsilon * m * Math.pow(sigma_ij, m) / Math.pow(r, m + 1);

        double force_mag = f_repulsive - f_attraction;

        double fx = force_mag * dx;
        double fy = force_mag * dy;

        ax[i] += fx;
        ay[i] += fy;
        ax[j] -= fx;
        ay[j] -= fy;

        totalPotentialEnergyAccumulator += epsilon * (Math.pow(sigma_ij / r, n) - Math.pow(sigma_ij / r, m));

        virialAccumulator += dx * fx + dy * fy;
      }
    }
  }

  private double pbcSeparation(double ds, double L) {
    if(ds>0) {
      while(ds>0.5*L) {
        ds -= L;
      }
    } else {
      while(ds<-0.5*L) {
        ds += L;
      }
    }
    return ds;
  }
  
  private double pbcPosition(double s, double L) {
    if(s>0) {
      while(s>L) {
        s -= L;
      }
    } else {
      while(s<0) {
        s += L;
      }
    }
    return s;
  }

  public void getRate(double[] state, double[] rate) {
    if(odeSolver.getRateCounter()==1) {
      computeAcceleration();
    }
    for(int i = 0;i<N;i++) {
      rate[4*i] = state[4*i+1];   
      rate[4*i+2] = state[4*i+3]; 
      rate[4*i+1] = ax[i];
      rate[4*i+3] = ay[i];
    }
    rate[4*N] = 1;
  }

  public double[] getState() {
    return state;
  }

  public void step(HistogramFrame xVelocityHistogram) {
    odeSolver.step();
    double totalKineticEnergy = 0;
    for(int i = 0; i < N; i++) {
      state[4*i+1] *= (1 - damping * dt); 
      state[4*i+3] *= (1 - damping * dt);

      totalKineticEnergy += (state[4*i+1]*state[4*i+1] + state[4*i+3]*state[4*i+3]);
      xVelocityHistogram.append(state[4*i+1]);
    }
    computeBoundingBox();
    totalKineticEnergy *= 0.5;
    steps++;
    totalKineticEnergyAccumulator += totalKineticEnergy;
    totalKineticEnergySquaredAccumulator += totalKineticEnergy * totalKineticEnergy;
    t += dt;
    double currentPE = totalPotentialEnergyAccumulator; 
    try (PrintWriter energyOut = new PrintWriter(new FileWriter("energy.txt", true))) {
      energyOut.printf("%d %.4f %.4f %.4f %.4f\n", steps, t, currentPE, width, height);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public void draw(DrawingPanel panel, Graphics g) {
    if(state==null) {
      return;
    }
    g.setColor(Color.red);
    for(int i = 0;i<N;i++) {
      int pxRadius = Math.abs(panel.xToPix(size[i])-panel.xToPix(0));
      int pyRadius = Math.abs(panel.yToPix(size[i])-panel.yToPix(0));
      int xpix = panel.xToPix(state[4*i])-pxRadius;
      int ypix = panel.yToPix(state[4*i+2])-pyRadius;
      g.drawOval(xpix, ypix, 2*pxRadius, 2*pyRadius);
    }
    g.setColor(Color.black);
    int xpix = panel.xToPix(0);
    int ypix = panel.yToPix(Ly);
    int lx = panel.xToPix(Lx)-panel.xToPix(0);
    int ly = panel.yToPix(0)-panel.yToPix(Ly);
    g.drawRect(xpix, ypix, lx, ly);
  }
}