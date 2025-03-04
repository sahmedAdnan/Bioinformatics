# common/utils.py
import os
import time
import shutil
from tqdm import tqdm
from colorama import Fore, Style, init

init(autoreset=True)

# Optional: Use pyfiglet for fancy banners.
try:
    from pyfiglet import Figlet
except ImportError:
    Figlet = None

contributors = ["Shakil", "Swantje",  "Ashrif"]

# Dynamically determine terminal width using shutil.
terminal_size = shutil.get_terminal_size(fallback=(80, 24))
DEFAULT_WIDTH = terminal_size.columns

class Colors:
    HEADER = Fore.CYAN + Style.BRIGHT
    OKBLUE = Fore.BLUE
    OKCYAN = Fore.CYAN
    OKGREEN = Fore.GREEN
    WARNING = Fore.YELLOW
    FAIL = Fore.RED
    MAGENTA = Fore.MAGENTA
    ENDC = Fore.RESET
    BOLD = Style.BRIGHT

def center_text(text, width=DEFAULT_WIDTH):
    return text.center(width)

def animate_banner(text, width=DEFAULT_WIDTH, delay=0.05):
    """
    Animate the banner by printing it line-by-line with a short delay.
    Uses the "big" font for a larger banner.
    """
    if Figlet:
        f = Figlet(font="digital")
        banner = f.renderText(text)
        lines = banner.splitlines()
        for line in lines:
            print(Colors.MAGENTA + line.center(width) + Colors.ENDC)
            time.sleep(delay)
    else:
        # Fallback: simple centered text
        print(Colors.HEADER + center_text("====== " + text + " ======", width) + Colors.ENDC)

def show_header(title="DNA Sequence Alignment ", width=DEFAULT_WIDTH):
    os.system("cls" if os.name == "nt" else "clear")
    animate_banner(title, width)
    print(Colors.HEADER + "=" * width)
    print(Colors.BOLD + center_text("👨‍💻Contributors: " + ", ".join(contributors), width))
    print(Colors.HEADER + "=" * width + Colors.ENDC)

def loading(message, duration=1.5, width=DEFAULT_WIDTH):
    print(Colors.WARNING + center_text(message, width))
    for _ in tqdm(range(50), bar_format="{l_bar}%s{bar}%s" % (Colors.OKGREEN, Colors.ENDC), leave=False):
        time.sleep(duration / 50)
    print(Colors.OKGREEN + center_text("[✔] Done!", width) + Colors.ENDC)

def validate_sequence(seq):
    valid_chars = {'A', 'C', 'G', 'T'}
    for char in seq:
        if char not in valid_chars:
            print(Colors.FAIL + f"✖ Error: Invalid character '{char}' found!" + Colors.ENDC)
            print(Colors.WARNING + "💡 Example of valid sequence: ACGTAGCTAGTC" + Colors.ENDC)
            return False
    return True

def get_int_input(prompt, error_message="✖ Please enter a valid integer."):
    while True:
        user_input = input(Colors.OKBLUE + prompt + Colors.ENDC)
        try:
            value = int(user_input)
            return value
        except ValueError:
            print(Colors.FAIL + error_message + Colors.ENDC)
