import pyautogui
import time

# Safety feature: moving the mouse to a corner will abort the script
pyautogui.FAILSAFE = True

silver_pens = True

# Coordinates for the three places to click (x, y)
silver_pen_points = [(130, 900), (200, 820), (200, 820)]

cog_points = [(100, 220), (180, 220), (260, 220),
              (100, 330), (180, 330), (260, 330),
              (100, 440), (180, 440), (260, 440),
              (100, 500), (180, 500), (260, 500),
              (100, 600), (180, 600), (260, 600),
              (280, 650)
              ]

# Delay to allow you to switch to the desired window
time.sleep(3)
if silver_pens:
    for _ in range(50000):
        for point in silver_pen_points:
            x, y = point
            pyautogui.moveTo(x, y, duration=0.001)  # Move to the specified coordinates
            pyautogui.click()  # Perform a left-click
else:
    for _ in range(4):
        for point in cog_points:
            x, y = point
            pyautogui.moveTo(x, y, duration=0.001)  # Move to the specified coordinates
            pyautogui.click()  # Perform a left-click

print("Script completed.")
