# Function to read the contacts from a given contact file and
# list of names and email addresses

def get_contacts(filename):
    names = []
    emails = []
    with open(filename, 'r') as contacts_file:
        for a_contact in contacts_file:
            names.append(a_contact.split()[0])
            emails.append(a_contact.split()[1])
    return names, emails

