"""Module embed_fonts.

Defines functionality useful for automatically embedding fonts in
svg files generated by Inkscape. The functionality is inspired by
the svgfontembed package (https://pypi.org/project/svgfontembed/).
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-05-04'
__license__ = 'GPLv3+'

from argparse import ArgumentParser
import base64
from collections import defaultdict
from contextlib import AbstractContextManager
import copy
from io import BytesIO
import re
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from fontTools.subset import Subsetter
from fontTools.ttLib import ttFont
from matplotlib.font_manager import findSystemFonts as sys_fonts
from parsel import Selector

FONT_FACE_RE = re.compile(r'@font-face\s*{[^}]*}', re.MULTILINE)
FONT_SRC_RE = re.compile(r'src:\s*url\(([^)]*)\)')
FONT_FAMILY_RE = re.compile(r'font-family:\s*([^;]*)')
FONT_WEIGHT_RE = re.compile(r'font-weight:\s*([^;]*)')
FONT_STYLE_RE = re.compile(r'font-style:\s*([^;\'"]*)')
_TAG_END = '>'
_MOD_FILE_NAME = '_embedded'


def load_system_fonts():
    """Return a list of FontFile(s) with all the known system fonts."""
    known_fonts = [Path(f) for f in sys_fonts()]
    fonts = [FontFile(f, ttFont.TTFont(f, fontNumber=0), 0)
             for f in known_fonts]
    # Each font file can contain multiple families. Collect those too.
    for font_file in fonts[:]:
        n_families = font_file.n_families
        path = font_file.path
        for i in range(1, n_families):
            fonts.append(FontFile(path, ttFont.TTFont(path, fontNumber=i), i))
    return fonts


@dataclass
class FontFace:
    """A class to represent a @font-face SVG definition.

    Example definition input (from the SVG):
        @font-face {
            font-family: "Virgil";
            src: url("https://somesite.com/Font.woff2");
            font-weight: bold;
            font-style: italic;
        }
    """

    definition: str
    # font_file_name: str | None = None

    @cached_property
    def family(self):
        """Return the font-family field of this font face."""
        font_name = FONT_FAMILY_RE.search(self.definition)
        if font_name:
            return font_name.group(1).strip('"\' ')
        return None

    @cached_property
    def src_url(self):
        """Return the contents of the url entry of the src field."""
        src = FONT_SRC_RE.search(self.definition)
        if src:
            return src.group(1).strip('"\' ')
        return None

    @cached_property
    def style(self):
        """Return the font-style field of this font face."""
        font_style = FONT_STYLE_RE.search(self.definition)
        if font_style:
            return font_style.group(1)
        return 'normal'

    @cached_property
    def weight(self):
        """Return the font-weight field of this font face."""
        font_weight = FONT_WEIGHT_RE.search(self.definition)
        if font_weight:
            return font_weight.group(1)
        return 'normal'

    @classmethod
    def from_attributes(cls, family, weight=None, style=None):
        """Return a FontFace from font attributes."""
        if not family.startswith(('"', '\'')):
            family = f'"{family}"'
        face_def = f'''  @font-face{{
      font-family: {family};
      src: url();
      font-weight: {weight or 'normal'};
      font-style: {style or 'normal'};
    }}
'''
        return cls(face_def)

    @classmethod
    def from_family_info(cls, info):
        """Return a FontFace from a FontFamilyInfo."""
        return cls.from_attributes(info.family, info.weight, info.style)

    @classmethod
    def from_svg(cls, svg_contents):
        """Return a FontFace(s) from @font-face definitions in `svg_contents`.

        Parameters
        ----------
        svg_contents : str
            The text read from an SVG file, from which @font-face
            definitions should be extracted.

        Returns
        -------
        tuple of FontFace
        """
        return tuple(cls(definition)
                     for definition in FONT_FACE_RE.findall(svg_contents))

    def with_woff_src(self, woff64):
        """Return a new FontFace with its src replaced by an encoded woff."""
        cls = type(self)
        src = f"src: url('data:font/woff2;base64,{woff64}') format('woff2');"
        result = re.sub(r'src:\s*url\(([^)]*)\)\s*(format(\'woff2\'))?\s*;',
                        src, self.definition)
        return cls(result)


@dataclass
class FontFile:
    """Container for one font face (including style) read from file.

    Attributes
    ----------
    path : Path
        Path to the file system from which the font was loaded.
    ttfont : fontTools.ttFont
        The font, loaded from file.
    font_num : int
        Which of the font faces contained in path this FontFile
        corresponds to. It only makes sense for font files that
        contain multiple font faces in the same file.
    """

    path: Path
    ttfont : ttFont
    font_num : int = 0

    @property
    def name(self):
        """Return the full font name of this font."""
        return self.ttfont['name'].getDebugName(4)

    @property
    def n_families(self):
        """Return the number of font families/faces in self.ttfont."""
        try:
            return self.ttfont.reader.numFonts
        except AttributeError:
            return 1

    def to_encoded_woff(self, subset=None):
        """Return the encoded glyphs, in woff2 format, as a string."""
        font = copy.deepcopy(self.ttfont)
        if subset is not None:
            subsetter = Subsetter()
            subsetter.populate(text=subset)
            subsetter.subset(font)
        font.flavor = 'woff2'
        bytestream = BytesIO()
        with BytesIO() as bytestream:
            font.save(bytestream)
            bytestream.seek(0)
            return base64.b64encode(bytestream.read()).decode('utf-8')


@dataclass
class FontFamilyInfo:
    """Container for a font family-weight-style definition from SVG <text>.

    Notice that this conceptually differs from FontFace in that the
    latter concerns a @font-face entry in a top-level <style> block,
    whereas FontFamilyInfo are collected from <text> and <tspan>
    blocks.
    """

    font_style_def: str

    @property
    def family(self):
        """Return the font-family field of this FontFamilyInfo."""
        match = FONT_FAMILY_RE.search(self.font_style_def)
        if match:
            return match.group(1)
        return None

    @family.setter
    def family(self, new_family):
        """Set a new font-family field for this FontFamilyInfo."""
        old_family = self.family
        if old_family is None:
            self.font_style_def += f'font-family:{new_family};'
            return
        self.font_style_def = FONT_FAMILY_RE.sub(f'font-family:{new_family}',
                                                 self.font_style_def)

    @property
    def full_name(self):
        """Return a reasonable full-font name from this FontFamilyInfo."""
        attrs = (getattr(self, attr) for attr in ('family', 'weight', 'style'))
        name = ' '.join(attr for attr in attrs if attr not in ('normal', None))
        return name.replace('\"', '').replace('\'', '')

    @property
    def style(self):
        """Return the font-style field of this FontFamilyInfo."""
        match = FONT_STYLE_RE.search(self.font_style_def)
        if match:
            return match.group(1)
        return None

    @style.setter
    def style(self, new_style):
        """Set a new font-style field for this FontFamilyInfo."""
        old_style = self.style
        if old_style is None:
            self.font_style_def += f'font-style:{new_style};'
            return
        self.font_style_def = FONT_STYLE_RE.sub(f'font-style:{new_style}',
                                                self.font_style_def)

    @property
    def weight(self):
        """Return the font-weight field of this FontFamilyInfo."""
        match = FONT_WEIGHT_RE.search(self.font_style_def)
        if match:
            return match.group(1)
        return None

    @weight.setter
    def weight(self, new_weight):
        """Set a new font-weight field for this FontFamilyInfo."""
        old_weight = self.weight
        if old_weight is None:
            self.font_style_def += f'font-weight:{new_weight};'
            return
        self.font_style_def = FONT_WEIGHT_RE.sub(f'font-weight:{new_weight}',
                                                 self.font_style_def)

    def update_from_other(self, other):
        """Replace fields in self with those of other that are defined."""
        for attr in ('family', 'weight', 'style'):
            other_attr = getattr(other, attr)
            if other_attr is not None:
                setattr(self, attr, other_attr)

    def __repr__(self):
        """Return a string representation of this FontFamilyInfo."""
        txt = ';'.join(f'font-{attr}:{getattr(self, attr)}'
                       for attr in ('family', 'weight', 'style')
                       if getattr(self, attr) is not None)
        return f'FontFamilyInfo("{txt};")'

    def __eq__(self, other):
        """Return whether this FontFamilyInfo is the same as other."""
        try:
            same_family = other.family == self.family
        except AttributeError:
            return NotImplemented
        return (same_family
                and (other.weight or 'normal') == (self.weight or 'normal')
                and (other.style or 'normal') == (self.style or 'normal'))

    def __hash__(self):
        """Return a hash from the attributes of this FontFamilyInfo."""
        optionals = self.weight or 'normal', self.style or 'normal'
        return hash((self.family, *optionals))


class FontEmbedder(AbstractContextManager):
    """A context manager that embeds fonts into an SVG file."""

    def __init__(self, path):
        """Initialize an instance from an SVG path."""
        self.path = Path(path)
        self.svg = ''

        # _families_and_characters: {FontFamilyInfo: set(characters)}
        self._families_and_characters = defaultdict(set)
        self._fonts = {}
        self._font_defs = []
        self._entered = False
        # pylint: disable-next=magic-value-comparison
        if self.path.suffix != '.svg':
            raise ValueError(f'{path} is not an SVG file')

    def __enter__(self):
        """Enter a context."""
        self._entered = True
        self._families_and_characters.clear()
        self._fonts.clear()
        self._font_defs.clear()
        self.svg = self.path.read_text(encoding='utf-8')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit context without suppressing any exception."""
        self._entered = False
        return super().__exit__(exc_type, exc_value, traceback)

    def embed_fonts(self, available_fonts):
        """Modify self.svg by embedding fonts."""
        if not self._entered:
            raise RuntimeError('Use only in a with statement')
        self._add_svg_style_tag()
        self._find_font_families()
        if not self._families_and_characters:
            print('Nothing to do: Found no text.')
            return
        self._collect_known_font_files(available_fonts)
        self._font_defs = FontFace.from_svg(self.svg)
        self._embed_fonts_impl()

    def save(self, path=None, rename=True, skip_if_useless=True):
        """Save modified SVG to a file path."""
        if (skip_if_useless
            and not self._families_and_characters
                and not self._font_defs):
            return
        if path is None and self.path.stem.endswith(_MOD_FILE_NAME):
            path = self.path
        elif path is None and rename:
            path = self.path.with_name(self.path.stem + _MOD_FILE_NAME
                                       + self.path.suffix)
        elif path is None:
            path = self.path
        with path.open('w', encoding='utf-8') as embedded:
            embedded.write(self.svg)

    def _add_svg_style_tag(self):
        """Add an empty <style></style> tag to self.svg` if not present yet."""
        has_style = Selector(self.svg).xpath('//svg/style').getall()
        if has_style:  # No need to add any
            return

        # Place an empty <style> field as the first child of <svg>
        svg_contents_lines = []
        found_svg_field = False
        for line in self.svg.splitlines():
            if line.strip().startswith('<svg'):
                found_svg_field = True
            if found_svg_field and _TAG_END in line:
                svg_end, *rest = line.split(_TAG_END, maxsplit=1)
                svg_contents_lines.append(svg_end + _TAG_END)
                svg_contents_lines.extend((
                    '  <style>',
                    '  </style>',
                    *rest,
                    ))
                found_svg_field = False
                continue
            svg_contents_lines.append(line)
        svg_mod = '\n'.join(svg_contents_lines)
        if not svg_mod.endswith('\n'):
            svg_mod += '\n'
        self.svg = svg_mod

    def _collect_known_font_files(self, known_fonts):
        """Select fonts that are suitable for any of `families`.

        Parameters
        ----------
        known_fonts : Iterable of FontFile
            The known fonts from which a selection should be extracted,
            containing only fonts that are pertinent to `families`.

        Returns
        -------
        font_files : dict
            Keys are FontFamilyInfo taken from `families`. Values are
            FontFile(s), taken from `known_fonts`, corresponding to
            the font file suitable for each family.
        """
        font_names = {f.name.lower(): f for f in known_fonts}
        for info in self._families_and_characters:
            search_name = info.full_name.lower()
            try:
                self._fonts[info] = next(f for fname, f in font_names.items()
                                         if fname == search_name)
            except StopIteration:
                pass
            else:
                continue
            # pylint: disable-next=magic-value-comparison
            if info.style == 'italic':
                # See if there's an oblique form of the same font
                search_name = search_name.replace('italic', 'oblique')
                try:
                    self._fonts[info] = next(
                        f for fname, f in font_names.items()
                        if fname == search_name
                        )
                except StopIteration:
                    pass
                else:
                    continue
            print(f'No font file available for {info}')

    def _embed_fonts_impl(self):
        """Actually edit self.svg with embedded fonts."""
        font_defs = self._font_defs
        for font_info, font_file in self._fonts.items():
            replace = True
            try:
                font_face = next(d for d in font_defs if font_info == d)
            except StopIteration:
                # Not found. Make a new one.
                replace = False
                font_face = FontFace.from_family_info(font_info)
            characters = ''.join(self._families_and_characters[font_info])
            woff_src = font_file.to_encoded_woff(characters)
            woff_def = font_face.with_woff_src(woff_src).definition
            if replace:
                self.svg = self.svg.replace(font_face.definition, woff_def)
            else:
                self.svg = self.svg.replace('</style>',
                                            f'{woff_def}  </style>')

    def _find_font_families(self):
        """Collect font families and characters used in self.svg."""
        # Font families may be specified in different spots:
        # in the @style of the top level <text> block or in
        # one of the @style(s) of the inner <tspan> blocks
        text_blocks = Selector(self.svg).xpath('.//text')
        for text_block in text_blocks:
            block_style = text_block.xpath('@style').getall()
            assert len(block_style) == 1
            block_info = FontFamilyInfo(block_style[0])
            self._walk_tspan_tree(text_block, block_info)

    def _walk_tspan_tree(self, node, info):
        """Update self._families_and_characters traversing `node` recursively.

        Parameters
        ----------
        node : parsel.Selector
            The parent node to start from. Commonly a <text> or <tspan>.
        info : FontFamilyInfo
            Font-face-style information already available in `node`.

        Returns
        -------
        None.
        """
        # The style information is structured in a hierarchical
        # manner: inner <tspan> blocks recursively inherit
        # unspecified font definitions from their parents
        for tspan in node.xpath('tspan'):
            tspan_info = copy.deepcopy(info)  # inherit from parent
            tspan_style = tspan.xpath('@style').getall()
            assert len(tspan_style) <= 1
            if tspan_style:
                this_info = FontFamilyInfo(tspan_style[0])
                tspan_info.update_from_other(this_info)
            characters = ''.join(tspan.xpath('text()').getall())
            if characters:
                self._families_and_characters[tspan_info].update(characters)
            self._walk_tspan_tree(tspan, tspan_info)


def embed_fonts_in_files(*files):
    """Embed fonts in a bunch of SVG files."""
    available_fonts = load_system_fonts()
    for file in files:
        print(f'Embedding fonts in {file}')
        with FontEmbedder(file) as svg:
            svg.embed_fonts(available_fonts)
            svg.save()


def main():
    """Embed fonts as per CLI arguments."""
    files = parse_args()
    embed_fonts_in_files(*files)


def parse_args():
    """Return an iterable of files to process from CLI arguments."""
    parser = ArgumentParser()
    parser.add_argument('--file', '-f')
    args = parser.parse_args()
    if args.file:
        return (args.file,)
    return Path(__file__).resolve().parents[2].glob('**/*.svg')


if __name__ == '__main__':
    main()
