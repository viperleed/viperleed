"""
======================================
  ViPErLEED Graphical User Interface
======================================
   *** Module guilib.widgetslib ***

Library of functions that are common to several Qt objects

Created: 2020-01-12
Author: Michele Riva
"""

import inspect
import re

import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed import guilib as gl
# import guilib as gl


################################################################################
#                                   FUNCTIONS                                  #
################################################################################


def editStyleSheet(qwidget, new_entries):
    """
    Append or replace new_entries in the styleSheet of qwidget
    
    Parameters
    ----------
    qwidget: QWidget
             QWidget whose style sheet will be edited
    new_entries: str
                 list of properties and values that are to be inserted in the
                 style sheet.
                 Format: "property1: value1; property2: value2; ..."
    """
    if not isinstance(qwidget, qtw.QWidget):
        raise ValueError("Argument 0 of editStyleSheet must be a QWidget "
                         f"subclass. Found {type(qwidget)} instead")

    if not isinstance(new_entries, str):
        raise TypeError("Argument 1 of editStyleSheet must be a string. "
                        f"Found {type(new_entries)} instead")

    entry_re = r"""
        (?P<property>[-a-z ]+)              # property:
                                            #    lowercase letters, '-' and ' '  
        :                                   # colon
        (?P<value>[\w\%\#\(\),:.-_\\\/ ]+)  # value:
                                            #    alphanumeric, %, #,
                                            #    parentheses, ' ' and
                                            #    .,:-_ \/
        """
    checkFormat = re.match(r"".join((r'^(', entry_re, r';)+$')),
                           new_entries, re.VERBOSE)
    if checkFormat is None:
        raise ValueError("Argument 1 of editStyleSheet must be a string "
                         "of the form 'property1: value1; "
                         "property2: value2; ...'")

    if len(qwidget.styleSheet()) == 0: # no style sheet present
        qwidget.setStyleSheet(new_entries)
        return

    # prepare a compiled re that will be used on all entries
    entry_format = re.compile(entry_re, re.VERBOSE)

    new_entries = (entry
                   for entry in re.split(';', new_entries)
                   if len(entry) > 0)

    for entry in new_entries:
        e_match = entry_format.match(entry)
        property = e_match.group('property')
        value = e_match.group('value')
        # check whether the entry is already present
        if property in qwidget.styleSheet():
            # update what's after "property:" till ";" (included) with value
            start = qwidget.styleSheet().index(property) + len(property) + 1
            # +1 to skip the ':'
            end = qwidget.styleSheet().index(';', start) + 1
            old = qwidget.styleSheet()
            qwidget.setStyleSheet("".join((old[:start], value, old[end:])))
        else: # otherwise simply append the new property
            qwidget.setStyleSheet(qwidget.styleSheet() + entry)


def drawText(painter, text, transform=None, combine=False):
    # reimplementation of text painting function that does not give awful
    # output. type(painter)==QPainter, type(p)==QPointF, type(text)==QString

    rawFont = qtg.QRawFont.fromFont(painter.font())
    indexes = rawFont.glyphIndexesForString(text)
    
    painter.save();
    paths = [rawFont.pathForGlyph(index) for index in indexes]
    advances = rawFont.advancesForGlyphIndexes(indexes,
                                               qtg.QRawFont.UseDesignMetrics
                                               | qtg.QRawFont.KernedAdvances)
    if transform is not None:
        painter.setWorldTransform(transform, combine=combine)
    for (path,advance) in zip(paths,advances):
        painter.fillPath(path, painter.pen().brush())
        painter.translate(advance)
    painter.restore()

def get_all_children_widgets(parent, exclude_parent=False, recursive=True):
    """
    This is an extension of the QObject.children() method that finds
    all the children QWidgets of a QObject.

    Parameters
    ----------
    parent: QWidget or QLayout (or subclasses)
    exclude_parent: bool, default = False
                    If True parent is not contained in the returned set
    recursive: bool, default = True
               If False, only first-generation widgets are returned 
               If True, all-generations widgets are included

    Returns
    -------
    list of QWidget
        - if parent is a QWidget:
          All QWidgets that are children of any depth of parent
        - if parent is a QLayout:
          All the QWidgets managed by parent and all their children of any depth
    """
    if not isinstance(parent, (qtw.QWidget, qtw.QLayout)):
        return set()
    
    children = set(parent.children())
    
    # find all widgets that are directly children of parent
    childrenWidgs = {child for child in [*children, parent]
                      if isinstance(child, qtw.QWidget)}
    # find all layouts that are children of parent
    childrenLays = {child for child in [*children, parent]
                      if isinstance(child, qtw.QLayout)}

    # add the widgets that are managed by the layouts in the list of children
    for lay in childrenLays:
        to_add = (lay.itemAt(idx).widget()
                  for idx in range(lay.count())
                  if isinstance(lay.itemAt(idx), qtw.QWidgetItem))
        childrenWidgs.update(widg
                             for widg in to_add
                             if isinstance(widg, qtw.QWidget))
    
    if recursive:
        # and run recursively to find all the nested children
        for child in childrenWidgs.copy():
            if child != parent:
                childrenWidgs.update(get_all_children_widgets(child))
    
    if exclude_parent:
        childrenWidgs.discard(parent)
    
    return childrenWidgs


################################################################################
#                                   CLASSES                                    #
################################################################################


class AllGUIFonts():  ## > Will handle in a different way!
    # May be worth removing the __init__ completely since all instances
    # have anyway the same 'value' for all fonts
    def __init__(self):
        self.buttonFont = qtg.QFont()
        self.smallTextFont = qtg.QFont()
        self.labelFont = qtg.QFont()
        self.smallButtonFont = qtg.QFont()
        self.largeTextFont = qtg.QFont()
        self.plotTitleFont = qtg.QFont()
        
        allfonts = [self.buttonFont,
                    self.smallTextFont,
                    self.labelFont,
                    self.smallButtonFont,
                    self.largeTextFont,
                    self.plotTitleFont]
        
        [x.setFamily("DejaVu Sans") for x in allfonts]
        
        self.buttonFont.setPointSize(10)
        self.smallTextFont.setPointSize(9)
        self.labelFont.setPointSize(10)
        self.smallButtonFont.setPointSize(8)
        self.largeTextFont.setPointSize(12)
        self.plotTitleFont.setPointSize(14)
        
        self.mathFont = qtg.QFont()
        self.mathFont.setFamily("CMU Serif")
        self.mathFont.setPointSize(15)
        self.mathFont.setStyleStrategy(qtg.QFont.PreferAntialias)
