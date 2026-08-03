/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package utils;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Print {
  private static final String highlight = "~~~~~~~~~~~~~~~~~~~~";

  public static final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

  public static void printInfo(String message) {
    printInfo(message, true);
  }

  public static void printInfo(String message, boolean lineFeed) {
    myPrint(message, "INFO", lineFeed);
  }

  public static void printInfoHighlight(String message) {
    message =  highlight + message + highlight;
    printInfo(message, true);
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
