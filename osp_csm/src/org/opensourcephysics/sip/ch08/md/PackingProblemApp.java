package org.opensourcephysics.sip.ch08.md;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.display.GUIUtils;


public class PackingProblemApp extends AbstractSimulation {
  PackingProblemParticles md = new PackingProblemParticles();
  PlotFrame pressureData = new PlotFrame("time", "PA/NkT", "Mean pressure");
  PlotFrame temperatureData = new PlotFrame("time", "temperature", "Mean temperature");
  HistogramFrame xVelocityHistogram = new HistogramFrame("vx", "H(vx)", "Velocity histogram");
  DisplayFrame display = new DisplayFrame("x", "y", "Lennard-Jones system");
  int T = 5000;

  public void initialize() {
    md.nx = control.getInt("nx");
    md.ny = control.getInt("ny");
    md.initialKineticEnergy = control.getDouble("initial kinetic energy per particle");
    md.Lx = control.getDouble("Lx");
    md.Ly = control.getDouble("Ly");
    md.initialConfiguration = control.getString("initial configuration");
    md.dt = control.getDouble("dt");
    md.initialize();
    display.addDrawable(md);
    display.setPreferredMinMax(0, md.Lx, 0, md.Ly);
    xVelocityHistogram.setBinWidth(2*md.initialKineticEnergy/md.N);
  }

  public void doStep() {
    if (md.steps >= T) {
      stop();
      return;
    }
    md.step(xVelocityHistogram);
    control.println("Step = " + md.steps);
    
    display.setMessage("Step = " + md.steps);

    try (PrintWriter out = new PrintWriter(new FileWriter("bounding_box_data.txt", true))) {
      if (md.steps == 1) {
        out.println("step time minX maxX minY maxY width height area");
      }
      out.printf("%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", 
          md.steps, md.t, md.minX, md.maxX, md.minY, md.maxY, md.width, md.height, md.area);
    } catch (IOException e) {
      e.printStackTrace();
    }
    pressureData.append(0, md.t, md.getMeanPressure());
    temperatureData.append(0, md.t, md.getMeanTemperature());
  }
  
  public void stop() {
    control.println("Density = "+decimalFormat.format(md.rho));
    control.println("Number of time steps = "+md.steps);
    control.println("Time step dt = "+decimalFormat.format(md.dt));
    control.println("<T>= "+decimalFormat.format(md.getMeanTemperature()));
    control.println("<E> = "+decimalFormat.format(md.getMeanEnergy()));
    control.println("Heat capacity = "+decimalFormat.format(md.getHeatCapacity()));
    control.println("<PA/NkT> = "+decimalFormat.format(md.getMeanPressure()));
  }
  
  public void startRunning() {
    md.dt = control.getDouble("dt");
    double Lx = control.getDouble("Lx");
    double Ly = control.getDouble("Ly");
    if((Lx!=md.Lx)||(Ly!=md.Ly)) {
      md.Lx = Lx;
      md.Ly = Ly;
      md.computeAcceleration();
      display.setPreferredMinMax(0, Lx, 0, Ly);
      resetData();
    }
  }

  public void reset() {
    control.setValue("nx", 6);
    control.setValue("ny", 6);
    control.setAdjustableValue("Lx", 20.0);
    control.setAdjustableValue("Ly", 15.0);
    control.setValue("initial kinetic energy per particle", 1.0);
    control.setAdjustableValue("dt", 0.01);
    control.setValue("initial configuration", "fixed");
    control.setValue("epsilon", 1.0);
    control.setValue("sigma", 1.0);
    control.setValue("n", 14);
    control.setValue("m", 4);
    control.setValue("damping", 0.1);   
    control.setValue("T", T);

    enableStepsPerDisplay(true);
    super.setStepsPerDisplay(10);
    display.setSquareAspect(true);
  }

  public void resetData() {
    md.resetAverages();
    GUIUtils.clearDrawingFrameData(false);
  }
  
  public static void main(String[] args) {
    SimulationControl control = SimulationControl.createApp(new PackingProblemApp());
    control.addButton("resetData", "Reset Data");
  }
}