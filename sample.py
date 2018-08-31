
class Sample():

    def __init__(self, id, variants = None, group = None):
        self.id = id
        if variants == None:
            self.variants = []
        else:
            self.variants = variants
        self.group = group