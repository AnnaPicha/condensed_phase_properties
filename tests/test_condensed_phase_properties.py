"""Test module for condensed_phase_properties."""

from condensed_phase_properties import __author__, __email__, __version__


def test_project_info():
    """Test __author__ value."""
    assert __author__ == "Anna Picha"
    assert __email__ == "anna.picha@gmail.com"
    assert __version__ == "0.0.0"
