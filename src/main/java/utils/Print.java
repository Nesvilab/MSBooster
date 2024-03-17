/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package utils;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Print {

  public static final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

  public static void printInfo(String message) {
    printInfo(message, true);
  }

  public static void printInfo(String message, boolean lineFeed) {
    myPrint(message, "INFO", lineFeed);
  }

  public static void printError(String message) {
    printError(message, true);
  }

  public static void printError(String message, boolean lineFeed) {
    myPrint(message, "ERROR", lineFeed);
  }

  private static void myPrint(String message, String level, boolean lineFeed) {
    if (lineFeed) {
      System.out.println(dateFormat.format(new Date(System.currentTimeMillis())) + " [" + level + "] - " + message);
    } else {
      System.out.print(dateFormat.format(new Date(System.currentTimeMillis())) + " [" + level + "] - " + message);
    }
  }

}
