""" from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from read_contacts import get_contacts
from read_template import read_template


def setup_smtp():
    import smtplib

    # set up the SMTP server

    server = smtplib.SMTP(host='smtp.cern.ch', port=587)
    server.starttls()
    server.login('mapse.b@cern.ch', 'MatterGravity45')
    return server


names, emails = get_contacts('mycontacts.txt')  # read contacts
message_template = read_template('msg.txt')

# For each contact, send the email:
for name, email in zip(names, emails):
    msg = MIMEMultipart() # create a message

    # add in the actual person name to the message template
    message = message_template.substitute(PERSON_NAME=name.title())

    # setup  the parameters of the message
    msg['From']='mapse.b@cern.ch'
    msg['To']=email
    msg['Subject']="This is a test!"

    # add in the message body
    msg.attach(MIMEText(message, 'plain'))

    # send the message via the server set up earlier
    server = setup_smtp()
    server.mail(msg)    

    del msg """

import smtplib

sender = 'mapse.b@cern.ch'
receivers = ['kevin.amarilo@cern.ch']

message = """From: From Mapse <mapse.b@cern.ch>
To: To Kelvin <kevin.amarilo@cern.ch>
Subject: SMTP e-mail test

This is a test e-mail message.
"""

server = smtplib.SMTP(host='smtp.cern.ch', port=587)
server.starttls()
server.login('mapse.b@cern.ch', 'MatterGravity45')    
server.sendmail(sender, receivers, message)         
print "Successfully sent email"

    
#['__doc__', '__ythoninit__', '__module__', '_get_socket', 'close', 'connect', 'data', 'debuglevel', 'default_port', 'docmd', 'does_esmtp', 'ehlo', 'ehlo_msg', 'ehlo_or_helo_if_needed', 'ehlo_resp', 'esmtp_features', 'expn', 'file', 'getreply', 'has_extn', 'helo', 'helo_resp', 'help', 'local_hostname', 'login', 'mail', 'noop', 'putcmd', 'quit', 'rcpt', 'rset', 'send', 'sendmail', 'set_debuglevel', 'sock', 'starttls', 'timeout', 'verify', 'vrfy']
