# for testing only, delete when done testing
import configparser
import csv
from time import localtime, strftime
from collections import defaultdict

def main():
    # config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
    # config.read('test_config.ini')
    # this_dict = config._sections['energy_calibration']
    this_dict = defaultdict(list)
    this_dict['nominal_energy'] = [1, 2, 3, 4, 5]
    name = 'test'
    clock = strftime("_%Y-%m-%d_%H-%M-%S", localtime())
    csv_name = name + clock + ".csv"
    with open(csv_name, 'w', encoding='UTF8', newline='') as file_name:
        writer = csv.writer(file_name)
        for key, value in this_dict.items():
            writer.writerow([key, value])

if __name__ == '__main__':
    main()
