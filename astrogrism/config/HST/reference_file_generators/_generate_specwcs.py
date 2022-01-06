from ast import literal_eval as eval
import datetime
import re
import numpy as np

import asdf
from astropy.modeling.models import Polynomial1D
from astropy import units as u

from astrogrism.HST.dispersion_models import DISPXY_Model, DISPXY_Extension
from astrogrism.config.HST.reference_file_generators._common import common_reference_file_keywords
from astrogrism.config.HST.reference_file_generators.wcs_ref_model import WFC3GrismModel


def dict_from_file(filename):
    """Read in a file and return a named tuple of the key value pairs.

    This is a generic read for a text file with the line following format:

    keyword<token>value

    Where keyword should start with a character, not a number
    Non-alphabetic starting characters are ignored
    <token> can be space or comma

    Parameters
    ----------
    filename : str
        Name of the file to interpret

    Examples
    --------
    dict_from_file('NIRCAM_C.conf')

    Returns
    -------
    dictionary of deciphered keys and values

    """
    # starts with a letter
    letters = re.compile("(^[a-zA-Z])")  # noqa: W605
    numbers = re.compile("(^(?:[+\-])?(?:\d*)(?:\.)?(?:\d*)?(?:[eE][+\-]?\d*$)?)")  # noqa: W605
    # is a blank line
    empty = re.compile("(^\s*$)")  # noqa: W605

    print("\nReading {0:s}  ...".format(filename))
    with open(filename, 'r') as fh:
        lines = fh.readlines()
    content = dict()
    for line in lines:
        value = None
        key = None
        if not empty.match(line):
            if letters.match(line):
                pair = re.split('\s+|(?<!\d)[,](?!\d)', line.strip())  # noqa: W605
                if len(pair) == 2:  # key and value exist
                    key = pair[0]  # first item is the key
                    val = pair[1]  # second item is the value
                    if letters.match(val):
                        value = val
                    if numbers.fullmatch(val):
                        value = eval(val)
                if len(pair) == 3:  # key min max exist
                    key = pair[0]
                    val1, val2 = pair[1:]
                    if numbers.fullmatch(val1) and numbers.fullmatch(val2):
                        value = (eval(val1), eval(val2))
                    else:
                        raise ValueError("Min/max values expected for {0}"
                                         .format(key))
                elif len(pair) > 3:  # Key and many values exist (DISPX and DISPY)
                    key = pair[0]
                    value = []
                    for i in range(1, len(pair)):
                        value.append(eval(pair[i]))
                    value = tuple(value)

        # ignore the filter file pointings and the sensitivity files these are
        # used for simulation
        if key and (value is not None):
            if (("FILTER" not in key) and ("SENSITIVITY" not in key)):
                content[key] = value

    return content


def split_order_info(keydict):
    """
    Designed to take as input the dictionary created by dict_from_file and for
    nircam, split out and accumulate the keys for each beam/order.
    The keys must have the beam in their string, the spurious beam designation
    is removed from the returned dictionary. Keywords with the same first name
    in the underscore separated string followed by a number are assumed to be
    ranges


    Parameters
    ----------
    keydict : dictionary
        Dictionary of key value pairs

    Returns
    -------
    dictionary of beams, where each beam has a dictionary of key-value pairs
    Any key pairs which are not associated with a beam get a separate entry
    """

    if not isinstance(keydict, dict):
        raise ValueError("Expected an input dictionary")

    # has beam name fits token
    token = re.compile('^[a-zA-Z]*_(?:[+\-]){0,1}[a-zA-Z0-9]{0,1}_*')  # noqa: W605
    rangekey = re.compile('^[a-zA-Z]*_[0-1]{0,1}[0-9]{1,1}$')  # noqa: W605
    rdict = dict()  # return dictionary
    beams = list()

    # prefetch number of Beams, beam is the second string
    for key in keydict:
        # Separate WEDGE keys that would otherwise match beam regex
        if key[0:5] == "WEDGE":
            if "WEDGE" not in beams:
                beams.append("WEDGE")
        elif token.match(key):
            b = key.split("_")[1].upper()
            if b not in beams:
                beams.append(b)
    for b in beams:
        rdict[b] = dict()

    #  assumes that keys are sep with underscore and beam is in second section
    for key in keydict:
        # Again, reject WEDGE keys
        if key[0:5] == "WEDGE":
            b, fname = key.split("_")
            rdict[b][fname] = keydict[key]
        elif token.match(key):
            b = key.split("_")[1].upper()
            newkey = key.replace("_{}_".format(b), "_")
            rdict[b][newkey] = keydict[key]

    # look for range variables to make them into tuples
    for b, d in rdict.items():
        keys = d.keys()
        rkeys = []
        odict = {}
        for k in keys:
            if rangekey.match(k):
                rkeys.append(k)
        for k in rkeys:
            mlist = [m for m in rkeys if k.split("_")[0] in m]
            root = mlist[0].split("_")[0]
            if root not in odict:
                temp_list = []
                for mk in mlist:
                    temp_list.append(d[mk])
                odict[root] = tuple(temp_list)
        # combine the dictionaries and remove the old keys
        d.update(odict)
        for k in rkeys:
            del d[k]

    return rdict


def create_grism_specwcs(conffile="",
                         pupil=None,
                         direct_filter=None,
                         author="STScI",
                         history="",
                         outname=None):
    """
    Note: This code is shamelessly stolen from the jwreftools package
    (see https://github.com/spacetelescope/jwreftools/) and adapted for use
    on HST GRISMCONF files. The docstrings and comments have not yet been
    updated accordingly.

    Create an asdf reference file to hold grism configuration information. No
    sensitivity information is included

    Note: The orders are named alphabetically, i.e. Order A, Order B
    There are also sensativity fits files which are tables of wavelength,
    sensativity, and error. These are specified in the conffile but will
    not be read in and saved in the output reference file for now.
    It's possible they may be included in the future, either here or as
    a separate reference files. Their use here would be to help define the
    min and max wavelengths which set the extent of the dispersed trace on
    the grism image. Convolving the sensitiviy file with the filter throughput
    allows one to calculate the wavelength of minimum throughput which defines
    the edges of the trace.

    direct_filter is not specified because it assumes that the wedge
    information (wx,wy) is included in the conf file in one of the key-value
    pairs, where the key includes the beam designation

     this reference file also contains the polynomial model which is
     appropriate for the coefficients which are listed.
     wavelength = DISPL(order,x0,y0,t)
     dx = DISPX(order,x0,y0,t)
     dy = DISPY(order,x0,y0,t)

     t = INVDISPX(order,x0,y0,dx)
     t = INVDISPY(order,x0,y0,dy)
     t = INVDISL(order,x0,y0, wavelength)



    Parameters
    ----------
    conffile : str
        The text file with configuration information, formatted as aXe expects
    pupil : str
        Name of the grism the conffile corresponds to
        Taken from the conffile name if not specified
    module : str
        Name of the Nircam module
        Taken from the conffile name if not specified
    author : str
        The name of the author
    history : str
        A comment about the refrence file to be saved with the meta information
    outname : str
        Output name for the reference file

    Returns
    -------
    fasdf : asdf.AsdfFile(WFC3GrismModel)

    """

    # if pupil is none get from filename like NIRCAM_modB_R.conf
    if pupil is None:
        pupil = "GRISM" + conffile.split(".")[0][-1]
    print("Pupil is {}".format(pupil))

    if outname is None:
        outname = "WFC3_{}_specwcs.asdf".format(pupil)
    if not history:
        history = "Created from {0:s}".format(conffile)

    if pupil in ("G102", "G141"):
        channel = "IR"
        wave_units = u.AA
    elif pupil == "G280":
        channel = "UVIS"
        wave_units = u.micron
    elif pupil == "G800L":
        channel = "WFC"
        wave_units = u.AA
    else:
        raise NotImplementedError("G102, G141, G280, G800L Grisms supported, not " + str(pupil))

    ref_kw = common_reference_file_keywords(reftype="specwcs",
                                            title=f"HST {channel} Grism Parameters",
                                            description="{0:s} dispersion models".format(pupil),
                                            exp_type=f"WFC3_{channel}",
                                            author=author,
                                            model_type="WFC3GrismModel",
                                            fname=direct_filter,
                                            pupil=pupil,
                                            filename=outname,
                                            )

    # get all the key-value pairs from the input file
    conf = dict_from_file(conffile)
    beamdict = split_order_info(conf)

    # Get x and y offsets from the filter, if necessary
    if direct_filter is not None:
        try:
            wx = beamdict["WEDGE"][direct_filter][0]
            wy = beamdict["WEDGE"][direct_filter][0]
        except KeyError:
            raise KeyError(f"No WEDGE information for {direct_filter} in conf file")
    else:
        wx, wy = 0, 0

    # UVIS conf files don't have WEDGE info
    try:
        beamdict.pop("WEDGE")
    except KeyError:
        pass

    # beam = re.compile('^(?:[+\-]){0,1}[a-zA-Z0-9]{0,1}$')  # match beam only
    # read in the sensitivity tables to save their content
    # they currently have names like this: NIRCam.A.1st.sensitivity.fits
    # translated as inst.beam/order.param
    temp = dict()
    # find beam key
    etoken = re.compile("^[a-zA-Z]*_(?:[+\-]){1,1}[1,2]{1,1}")  # noqa: W605
    for b, bdict in beamdict.items():
        temp[b] = dict()

    # add the new beam information to beamdict and remove spurious beam info
    for k in temp:
        for kk in temp[k]:
            if etoken.match(kk):
                kk = kk.replace("_{}".format(k), "")
            beamdict[k][kk] = temp[k][kk]

    # for NIRCAM, the R and C grism coefficients contain zeros where
    # the dispersion is in the opposite direction. Meaning, the GRISMR,
    # which disperses along ROWS has coefficients of zero in the y models
    # and vice versa.
    #
    # There are separate reference files for each grism. Depending on the grism
    # dispersion direction you either want to use the dx from source center or
    # the dy from source center in the inverse dispersion relationship which is
    # used to calculate the t value needed to calculate the wavelength at that
    # pixel.
    # The model creation here takes all of this into account by looking at the
    # GRISM[R/C] the file is used for and creating a reference model with the
    # appropriate dispersion direction in use. This eliminates having to decide
    # which direction to calculatethe dispersion from given the input x,y
    # pixel in the dispersed image.
    orders = beamdict.keys()

    print(f"Orders: {orders}")

    # dispersion models valid per order and direction saved to reference file
    # Forward
    invdispl = []
    invdispx = []
    invdispy = []
    # Backward
    displ = []
    dispx = []
    dispy = []

    for order in orders:
        # convert the displ wavelengths from Angstroms (assumed default)
        # to microns for not IR grisms G800L and G280
        l_coeffs = np.array(beamdict[order]['DISPL'])
        if wave_units == u.micron:
            l_coeffs = (l_coeffs * u.AA).to_value(u.micron)

        # create polynomials using the coefficients of each order

        # This holds the wavelength lookup coeffs
        # This model is  INVDISPL for backward and returns t
        # This model should be DISPL for forward and returns wavelength
        if l_coeffs.shape == (2,):
            l0 = l_coeffs[0]
            l1 = l_coeffs[1]
            if l1 == 0:
                lmodel = Polynomial1D(1, c0=0, c1=0)
            else:
                lmodel = Polynomial1D(1, c0=-l0/l1, c1=1./l1)
            invdispl.append(lmodel)
            lmodel = Polynomial1D(1, c0=l0, c1=l1)
            displ.append(lmodel)
        else:
            # Can't invert higher order in t
            if len(l_coeffs.shape) > 1:
                try:
                    lmodel = DISPXY_Model(l_coeffs, 0, inv=True)
                except ValueError:
                    lmodel = None
            invdispl.append(lmodel)

            lmodel = DISPXY_Model(l_coeffs, 0)
            displ.append(lmodel)

        # This holds the x coefficients, for the R grism this model is the
        # the INVDISPX returning t, for the C grism this model is the DISPX
        e = beamdict[order]['DISPX']
        xmodel = DISPXY_Model(e, wx)
        dispx.append(xmodel)
        try:
            inv_xmodel = DISPXY_Model(e, wx, inv=True)
            invdispx.append(inv_xmodel)
        except Exception:
            # Interpolate x inverse
            invdispx.append(None)

        # This holds the y coefficients, for the C grism, this model is
        # the INVDISPY, returning t, for the R grism, this model is the DISPY
        e = beamdict[order]['DISPY']
        ymodel = DISPXY_Model(e, wy)
        dispy.append(ymodel)
        try:
            inv_ymodel = DISPXY_Model(e, wy, inv=True)
            invdispy.append(inv_ymodel)
        except Exception:
            invdispy.append(None)

    # We need to register the converter for the DISPXY_Model class with asdf
    asdf.get_config().add_extension(DISPXY_Extension())

    ref = WFC3GrismModel(channel=channel)
    ref.meta.update(ref_kw)
    # This reference file is good for NRC_WFSS and TSGRISM modes
    ref.meta.exposure.p_exptype = "NRC_WFSS|NRC_TSGRISM"
    ref.meta.input_units = u.micron
    ref.meta.output_units = u.micron
    ref.displ = displ
    ref.dispx = dispx
    ref.dispy = dispy
    ref.invdispx = invdispx
    ref.invdispy = invdispy
    ref.invdispl = invdispl
    # change the orders into translatable integers
    # so that we can look up the order with the proper index
    ref.order = [int(o) for o in beamdict]
    history = asdf.tags.core.HistoryEntry({'description': history,
                                           'time': datetime.datetime.utcnow()})
    software = asdf.tags.core.Software({'name': '_generate_specwcs.py',
                                        'author': author,
                                        'homepage': 'https://github.com/spacetelescope/astrogrism',
                                        'version': '1.0.0'})

    history['software'] = software
    ref.history = [history]
    ref.to_asdf(outname)
    ref.validate()
