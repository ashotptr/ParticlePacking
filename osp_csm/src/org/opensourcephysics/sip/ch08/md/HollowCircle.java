package org.opensourcephysics.sip.ch08.md;

import java.awt.Color;
import java.awt.Graphics;

import org.opensourcephysics.display.Circle;
import org.opensourcephysics.display.DrawingPanel;

public class HollowCircle extends Circle {
    public HollowCircle() {
        super();
        this.color = Color.BLACK; // Set default edge color
    }
  
    @Override
    public void draw(DrawingPanel panel, Graphics g) {
        int xpix = panel.xToPix(this.x) - this.pixRadius;
        int ypix = panel.yToPix(this.y) - this.pixRadius;
        g.setColor(this.color);
        // Draw outline instead of filled circle
        g.drawOval(xpix, ypix, 2 * this.pixRadius, 2 * this.pixRadius);
    }
  }
  