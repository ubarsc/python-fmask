:: This file is part of 'python-fmask' - a cloud masking module
:: Copyright (C) 2015  Neil Flood
::
:: This program is free software; you can redistribute it and/or
:: modify it under the terms of the GNU General Public License
:: as published by the Free Software Foundation; either version 2
:: of the License, or (at your option) any later version.
::
:: This program is distributed in the hope that it will be useful,
:: but WITHOUT ANY WARRANTY; without even the implied warranty of
:: MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
:: GNU General Public License for more details.
::
:: You should have received a copy of the GNU General Public License
:: along with this program; if not, write to the Free Software
:: Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

:: Run a Python script in the same dir as this, and same name, but .py extension
:: with the same command line args
@echo off
python %~dpn0.py %*
