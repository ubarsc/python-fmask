"""
Exceptions used within fmask
"""

# This file is part of 'python-fmask' - a cloud masking module
# Copyright (C) 2015  Neil Flood
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


class FmaskException(Exception):
    "An exception rasied by Fmask"


class FmaskParameterError(FmaskException):
    "Something is wrong with a parameter"


class FmaskFileError(FmaskException):
    "Data in file is incorrect"


class FmaskNotSupportedError(FmaskException):
    "Requested operation is not supported"


class FmaskInstallationError(FmaskException):
    "Problem with installation of Fmask or a required package"


class Sen2MetaError(Exception):
    "Error with Sentinel-2 metadata"
