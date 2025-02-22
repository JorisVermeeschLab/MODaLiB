from modalibpy.regions.region import Region


class Gene(Region):

    def __init__(self, name: str, *args) -> None:
        Region.__init__(self, *args)
        self.name: str = str(name)

    def __str__(self) -> str:
        return f'Gene({self.name})'

    def __repr__(self) -> str:
        return self.__str__()
