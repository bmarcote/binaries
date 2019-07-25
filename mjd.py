#!/usr/bin/env python
import datetime
from sys import argv
import calendar
import time
from numpy import fmod

def date2mjd(date):
  origin = datetime.datetime(1858,11,17)
  mjd = (date-origin).days + (date-origin).seconds/86400.0
  return mjd

def mjd2date(mjd):
  origin = datetime.datetime(1858,11,17)
  date = origin + datetime.timedelta(mjd)
  return date

def mjd2jd(date):
    return date + 2400000.5

def jd2mjd(date):
    return date - 2400000.5

def date2year(date):
    s = lambda date: calendar.timegm(date.timetuple())
    year = date.year
    startOfThisYear = datetime.datetime(year=year, month=1, day=1)
    startOfNextYear = datetime.datetime(year=year+1, month=1, day=1)
    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = 1.0*yearElapsed/yearDuration
    return date.year + fraction

def year2date(year):
    s = lambda date: calendar.timegm(date.timetuple())
    fraction = fmod(year, 1)
    startOfThisYear = datetime.datetime(year=int(year), month=1, day=1)
    startOfNextYear = datetime.datetime(year=int(year)+1, month=1, day=1)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    yearElapsed  = fraction*yearDuration
    date_struc = time.gmtime(s(startOfThisYear) + yearElapsed)
    date = datetime.datetime.utcfromtimestamp(calendar.timegm(date_struc))
    return date

def prt_date(d):
    mjd, date = datetomjd(d)
    print(date, date.ctime()[:7])
    print('mjd = %11.5f' %mjd)

if __name__ == '__main__':
    if (len(argv) != 7) and (len(argv) != 2):
        print(""" Error, introduce 6 o 1 parametros:
            Para obtener el MJD:      mjd.py d m y h m s
            Para obtener la fecha:    mjd.py mjd""")
        exit(1)
    try:
        day = int(argv[1])
        month = int(argv[2])
        year = int(argv[3])
        hour = int(argv[4])
        minut = int(argv[5])
        second= int(argv[6])
        date = datetime.datetime(year, month, day, hour, minut, second)    
        prt_date(date)
    except:
        mjd = float(argv[1])
        date = mjd2date(mjd)
        prt_date(date)

