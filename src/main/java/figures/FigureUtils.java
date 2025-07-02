package figures;

import java.awt.*;

public class FigureUtils {

    //suggested by ChatGPT
    static Color getColor(int index) {
        // Hue cycles from 0 to 1 (full spectrum)
        float hue = (index * 0.61803398875f) % 1; // Golden ratio for better spacing
        return Color.getHSBColor(hue, 0.7f, 0.95f);
    }

    static Color getColorWithAlpha(int index, float alpha) {
        float hue = (index * 0.61803398875f) % 1.0f;  // Golden ratio hue spacing
        float saturation = 0.7f;
        float brightness = 0.95f;

        int rgb = Color.HSBtoRGB(hue, saturation, brightness);
        int r = (rgb >> 16) & 0xFF;
        int g = (rgb >> 8) & 0xFF;
        int b = rgb & 0xFF;
        int a = Math.round(alpha * 255);  // alpha in range [0.0, 1.0]

        return new Color(r, g, b, a);
    }

}
