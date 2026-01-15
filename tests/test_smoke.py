from dig_eco import hello


def test_hello_smoke() -> None:
    assert hello("world") == "hello, world"
