from density import MainWindow


def test_warning_banner_shows_cleaned_messages_and_restores_layout(qt_app):
    window = MainWindow(lhv_data={})
    base_results_label_geometry = window.ui.results_label.geometry()
    base_results_list_geometry = window.ui.results_list.geometry()
    base_progress_geometry = window.ui.progressBar.geometry()

    window.set_thermo_warning_messages(["  first warning  ", "", "second"])

    assert not window.thermo_warning_label.isHidden()
    assert window.thermo_warning_label.text() == "first warning\nsecond"
    assert window.ui.results_label.geometry().y() > (
        base_results_label_geometry.y()
    )
    assert window.ui.results_list.geometry().height() <= (
        base_results_list_geometry.height()
    )

    window.set_thermo_warning_messages([])

    assert window.thermo_warning_label.isHidden()
    assert window.ui.results_label.geometry() == base_results_label_geometry
    assert window.ui.results_list.geometry() == base_results_list_geometry
    assert window.ui.progressBar.geometry() == base_progress_geometry

    window.close()
